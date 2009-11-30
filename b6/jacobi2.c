#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "timer.h"
#include <mpi.h>

int rank;

// Maximums-Norm
double max_norm(const double* const x, const double* const b, const double* const A, long n)
{
  double max = 0.0;
  long i,j;

  for (i=0; i<n; i++)
  {
    double res = b[i];
    for (j=0; j<n; j++)
      res -= A[i*n+j] * x[j];
    if (max < fabs(res))
      max  = fabs(res);
  }

  return max;
}

// Initialisierung
void init(double* x, double* A, long n)
{
  long i,j;

  for (i=0; i<n; i++)
  {
    x[i] = 0.0;
    for (j=0; j<n; j++)
      A[i*n+j] = 0.0;

    A[i*n+i] = 3.0;
  }
}

void initB(double* b, long n)
{
	int c;
	for(c=0;c<n;c++){
		b[c]=1.0;
	}
}

// checke einen Vektor
void output(double* x, long n)
{
  long i;
  for (i=0; i<n; i++)
    fprintf(stdout, "%.4e ", x[i]);
  fprintf(stdout, "\n");
}

// Jacobi-Schritt
void jacobi(double* xneu, const double* const xalt, const double* const b, const double* const A, long n)
{
  printf("R(%d):\n",rank);
  printf("xneu0(%f)\n",xneu[0]);
  printf("xneu1(%f)\n",xneu[1]);
  printf("xalt0(%f)\n",xalt[0]);
  printf("xalt1(%f)\n",xalt[1]);
  printf("b0(%f)\n",b[0]);
  printf("b1(%f)\n",b[1]);
  printf("A00(%f)\n",A[0]);
  printf("A01(%f)\n",A[1]);
  printf("n(%d)\n",n);

  int i=0, j=0;

  for (i=0; i<n; ++i)
  {
    xneu[i] = 0.0;

    // compute scalar product
    for (j=0; j<n; ++j)
      xneu[i] += A[i*n+j] * xalt[j];

    // subtract case i=j again
    xneu[i] -= A[i*n+i] * xalt[i];

    // finish jacobi-step
    xneu[i] *= -1.0;
    xneu[i] += b[i];
    if (A[i*n+i] == 0)
    {
      fprintf(stderr, "Division by 0\n");
      exit (1);
    }
    else
      xneu[i] /= A[i*n+i];
  }
}

void calcSendcounts(int n, int psize, int *scount_A, int *displs_A)
{
	int partition = n/psize;
	int c, offset_A=0,offset_b=0;
	for(c=0;c<psize;c++){
		scount_A[c]=partition*n;
		displs_A[c]=offset_A*n;
		offset_A+=partition;
	}
	if(n%psize!=0){
		scount_A[psize-1]+=(n%psize)*n;
	}
}

void calcHelperMatrix(int n, int rank, int psize, double** B)
{
	int partition=n/psize;
	if(n%psize!=0 && rank==psize-1)partition+=n%psize;
	*B = (double *) malloc(sizeof(double)*n*partition);
}

//*****************************************************************************
// main
//*****************************************************************************
int main(int argc, char **argv)
{
  const int maxIt = 10000; // max Iterationen des Verfahrens
  long it;                 // Schleifenzaehler
  long n;                  // Problemgroesse in einer Richtung
  int psize;
  int *sendcounts_A;
  int *senddispls_A;
//  int rank;
  double *x;               // Unbekannte
  double *y;               // alte Unbekannte
  double *b;               // rechte Seite
  double *A;               // Matrix
  double *B;               // Matrix
  struct timeval timer;
  const double eps = 1e-5; // Abbruchbedingung
  double diff;             // Maximumsnorm des Abstandes zweier Loesungen
  double t = 0.0;          // fuer die Zeitmessung
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&psize);
  // Lese Laenge des Unbekannten-Vektors x
  if (argc != 2)
  {
    fprintf(stderr, "Syntax: %s <n>\n", argv[0]);
    exit (1);
  }
  n = atoi(argv[1]);
  sendcounts_A = (int *) malloc(sizeof(int)*psize);
  senddispls_A = (int *) malloc(sizeof(int)*psize);
  calcSendcounts(n,psize,sendcounts_A,senddispls_A);
  calcHelperMatrix(n,rank,psize,&B);
  x = (double *) malloc(sizeof(double)*n);
  y = (double *) malloc(sizeof(double)*n);
  if(rank==0){
    A = (double *) malloc(sizeof(double)*n*n);
  }
  b = (double *) malloc(sizeof(double)*n);
  // Initialisiere A, b und x
  if(rank==0)init(x, A, n);
  initB(b,n);
  // Beginn Zeitmessung
  if(rank==0)reset_timer(&timer);

  // Jacobi-Schritte
  for (it=1; it<((long) (maxIt/2)+1); it++)
  {
    MPI_Scatterv(A,sendcounts_A,senddispls_A,MPI_DOUBLE,B,sendcounts_A[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(y,x,b,B,sendcounts_A[rank]);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(B,sendcounts_A[rank],MPI_DOUBLE,A,sendcounts_A,senddispls_A,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(A,sendcounts_A,senddispls_A,MPI_DOUBLE,B,sendcounts_A[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(x,y,b,B,sendcounts_A[rank]);
    MPI_Barrier(MPI_COMM_WORLD);
    //output(x, n);

    if(rank==0){
      // Konvergenz-Check und Ende Zeitmessung
      diff = max_norm(x,b,A,n);
      fprintf(stdout, "Norm(res): %.16e\n", diff);
      if (diff < eps)
      {
        t = get_timer(timer);
        fprintf(stdout, "Iter: %d, t: %f s, Norm(res): %.16e\n", (2*it), t, diff);
        break;
      }
    }
  }

  // release memory
  free(x);
  free(y);
  free(b);
  free(A);
  free(B);
  MPI_Finalize();
  return (0);
}
