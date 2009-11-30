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

// Jacobi-Schritt
void jacobi(double* xneu, const double* const xalt, const double* const b, const double* const A, long start, long width, long n)
{
  int i=0, j=0;

  for (i=0; i<width; ++i)
  {
    xneu[i] = 0.0;

    // compute scalar product
    for (j=0; j<n; ++j)
      xneu[i] += A[(i+start)*n+j] * xalt[j];

    // subtract case i=j again
    xneu[i] -= A[(i+start)*n+(i+start)] * xalt[i+start];

    // finish jacobi-step
    xneu[i] *= -1.0;
    xneu[i] += b[i+start];
    if (A[(i+start)*n+(i+start)] == 0)
    {
      fprintf(stderr, "Division by 0\n");
      exit (1);
    }
    else
      xneu[i] /= A[(i+start)*n+(i+start)];
  }
}

void calcSendCounts(int *sendcounts, int *senddispls, int n, int psize)
{
	int partition = ceil(1.0*n/psize);
	int c, offset;
	for(c=0;c<psize;c++){
		sendcounts[c] = partition;
		senddispls[c] = offset;
		offset += partition;
	}
	if(n%psize!=0)sendcounts[psize-1]+=n%psize;
}

//*****************************************************************************
// main
//*****************************************************************************
int main(int argc, char **argv)
{
  const int maxIt = 10000; // max Iterationen des Verfahrens
  long it;                 // Schleifenzaehler
  long n;                  // Problemgroesse in einer Richtung

  double *x;               // Unbekannte
  double *y;               // alte Unbekannte
  double *b;               // rechte Seite
  double *A;               // Matrix
  double *xx;
  double *yy;

  int *sendcounts;
  int *senddispls;
  int rank;
  int psize;

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

  // allokiere Speicher fuer A, x, y und b
  x = (double *) malloc(sizeof(double)*n);   // alter Vektor x
  y = (double *) malloc(sizeof(double)*n);   // neuer Vektor y
  b = (double *) malloc(sizeof(double)*n);   // rechte Seite b
  A = (double *) malloc(sizeof(double)*n*n); // Matrix A
  sendcounts = (int *) malloc(sizeof(int)*psize);   // alter Vektor x
  senddispls = (int *) malloc(sizeof(int)*psize);   // alter Vektor x

  // Initialisiere A, b und x
  init(x, b, A, n);

  calcSendCounts(sendcounts,senddispls,n,psize);
  xx = (double *) malloc(sizeof(double)*sendcounts[rank]);
  yy = (double *) malloc(sizeof(double)*sendcounts[rank]);

  // Beginn Zeitmessung
  struct timeval timer;
  if(rank==0)reset_timer(&timer);

  // Jacobi-Schritte
  for (it=1; it<((long) (maxIt/2)+1); it++)
  {
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(yy,x,b,A,rank*sendcounts[0],sendcounts[rank],n);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(yy,sendcounts[rank],MPI_DOUBLE,y,sendcounts,senddispls,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(y,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(xx,y,b,A,rank*sendcounts[0],sendcounts[rank],n);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(xx,sendcounts[rank],MPI_DOUBLE,x,sendcounts,senddispls,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

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
  free(xx);
  free(yy);
  free(b);
  free(A);
  MPI_Finalize();
  return (0);
}
