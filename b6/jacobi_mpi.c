#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "timer.h"
#include<mpi.h>

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
    b[i] = i;
    for (j=0; j<n; j++)
      A[i*n+j] = 2.0; // hier gibts nen Fehler

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
void jacobi(double* xneu, const double* const xalt, const double* const b, const double* const A, long n, int rank, int procs)
{
  long i=0, j=0;
  long x=(long)ceil(1.0*n/procs);

  for (i=rank*x; i<fmin((rank+1)*x,n); ++i)
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

bool tolerance_reached(double* x, double* b, double* A, long n){
  const double eps = 1e-5; // Abbruchbedingung
  double diff;

  diff = max_norm(x,b,A,n);
  if(false)
    fprintf(stdout, "Norm(res): %.16e\n", diff);

  return diff < eps;
}
void d(int i){
  printf("step %i\n",i);
}

//*****************************************************************************
// main
//*****************************************************************************
int main(int argc, char **argv)
{
  const int maxIt = 10000; // max Iterationen des Verfahrens
  long it;                 // Schleifenzaehler
  long n;                  // Problemgroesse in einer Richtung
  int i;
  long interval;

  double *x;               // Unbekannte
  double *y;               // alte Unbekannte
  double *b;               // rechte Seite
  double *A;               // Matrix

  double diff;             // Maximumsnorm des Abstandes zweier Loesungen
  double t = 0.0;          // fuer die Zeitmessung

  int rank, procs;
  struct timeval timer;
  
  MPI_Status status;
  MPI_Init(&argc, &argv);

  /* Get own rank and total number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // Lese Laenge des Unbekannten-Vektors x
  if (argc != 2)
  {
    fprintf(stderr, "Syntax: %s <n>\n", argv[0]);
    exit (1);
  }
  n = atol(argv[1]);

  // allokiere Speicher fuer A, x, y und b
  x = (double *) malloc(sizeof(double)*n);   // alter Vektor x
  y = (double *) malloc(sizeof(double)*n);   // neuer Vektor y
  b = (double *) malloc(sizeof(double)*n);   // rechte Seite b
  A = (double *) malloc(sizeof(double)*n*n); // Matrix A

  // Initialisiere A, b und x
  init(x, b, A, n);

  if(!rank){
    // Beginn Zeitmessung
    reset_timer(&timer);
  }
  it=0;

  while(!tolerance_reached(x,b,A,n) && it<((long) (maxIt))){

    jacobi(y,x,b,A,n,rank,procs);
    //output(x, n);

    interval=(long)ceil(1.0*n/procs);

    for(i=0;i<procs;i++)
      if((i+1)*interval>n){
	MPI_Bcast(&y[i*interval], n-i*interval, MPI_DOUBLE, i, MPI_COMM_WORLD );
	break;
      }
      else{
	MPI_Bcast(&y[i*interval], interval, MPI_DOUBLE, i, MPI_COMM_WORLD );
      }

    x=y;
    ++it;
  }

  if(!rank){
    t = get_timer(timer);
    fprintf(stdout, "Iter: %d, t: %f s, Norm(res): %.16e\n", it+1, t, max_norm(x,b,A,n));
  }

  MPI_Finalize();

  // release memory
  free(x);
  free(b);
  free(A);

  return (0);
}
