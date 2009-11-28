//******************************************************************************
// jacobi-seq.c
//
// Jacobi-Iteration mit (vollbesetzter) Matrix, sequentielle Variante
//******************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "timer.h"
#include <mpi.h>

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
void init(double* x, double* b, double* A, long n)
{
  long i,j;

  for (i=0; i<n; i++)
  {
    x[i] = 0.0;
    b[i] = 1.0;
    for (j=0; j<n; j++)
      A[i*n+j] = 0.0;

    A[i*n+i] = 3.0;
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
      fprintf(stderr, "Division by 0\n.");
      exit (1);
    }
    else
      xneu[i] /= A[i*n+i];
  }
}

void calcSendcounts(int n, int psize, int *scount_A, int *scount_b, int *displs_A, int *displs_b)
{
	int partition = n/psize;
	int c, offset_A=0,offset_b=0;
	for(c=0;c<psize;c++){
		scount_A[c]=partition*n;
		scount_b[c]=partition;
		displs_A[c]=offset_A*n;
		displs_b[c]=offset_b;
		offset_A+=partition;
		offset_b+=partition;
	}
	if(n%psize!=0){
		scount_A[psize-1]+=(n%psize)*n;
		scount_b[psize-1]+=n%psize;
	}
}

void calcHelperMatrix(int n, int rank, int psize, double** B, double** c)
{
	int partition=n/psize;
	if(n%psize!=0 && rank==psize-1)partition+=n%psize;
	*B = (double *) malloc(sizeof(double)*n*partition);
	*c = (double *) malloc(sizeof(double)*partition);
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
  int *sendcounts_b;
  int *senddispls_A;
  int *senddispls_b;
  int rank;
  double *x;               // Unbekannte
  double *y;               // alte Unbekannte
  double *b;               // rechte Seite
  double *c;               // rechte Seite
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

  x = (double *) malloc(sizeof(double)*n);
  y = (double *) malloc(sizeof(double)*n);
  c = (double *) malloc(sizeof(double)*n);
  sendcounts_A = (int *) malloc(sizeof(int)*psize);
  sendcounts_b = (int *) malloc(sizeof(int)*psize);
  senddispls_A = (int *) malloc(sizeof(int)*psize);
  senddispls_b = (int *) malloc(sizeof(int)*psize);
  calcSendcounts(n,psize,sendcounts_A,sendcounts_b,senddispls_A,senddispls_b);
  calcHelperMatrix(n,rank,psize,&B,&c);
  if(rank==0){
    A = (double *) malloc(sizeof(double)*n*n);
    b = (double *) malloc(sizeof(double)*n);
  }
  // Initialisiere A, b und x
  if(rank==0)init(x, b, A, n);
  // Beginn Zeitmessung
  if(rank==0)reset_timer(&timer);

  // Jacobi-Schritte
  for (it=1; it<((long) (maxIt/2)+1); it++)
  {
    MPI_Scatterv(A,sendcounts_A,senddispls_A,MPI_DOUBLE,B,sendcounts_A[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatterv(b,sendcounts_b,senddispls_b,MPI_DOUBLE,c,sendcounts_b[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(y,x,c,B,n);
    jacobi(x,y,c,B,n);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(B,sendcounts_A[rank],MPI_DOUBLE,A,sendcounts_A,senddispls_A,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gatherv(c,sendcounts_b[rank],MPI_DOUBLE,b,sendcounts_b,senddispls_b,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
  free(c);
  free(A);
  free(B);
  MPI_Finalize();
  return (0);
}
