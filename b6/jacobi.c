//******************************************************************************
// jacobi-seq.c
//
// Jacobi-Iteration mit (vollbesetzter) Matrix, sequentielle Variante
//******************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "timer.h"

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
      A[i*n+j] = 1.0;

    A[i*n+i] = n+1;
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

bool tolerance_reached(const double* x, const double* b, const double* A, const long n){
  const double eps = 1e-5; // Abbruchbedingung
  double diff;

  diff = max_norm(x,b,A,n);
  if(false)
    fprintf(stdout, "Norm(res): %.16e\n", diff);

  return diff < eps;
}

//*****************************************************************************
// main
//*****************************************************************************
int main(int argc, char **argv)
{
  const int maxIt = 10000; // max Iterationen des Verfahrens
  long it;                 // Schleifenzaehler
  long n;                  // Problemgroesse in einer Richtung
  long i;

  double *x;               // Unbekannte
  double *y;               // alte Unbekannte
  double *b;               // rechte Seite
  double *A;               // Matrix

  const double eps = 1e-5; // Abbruchbedingung
  double diff;             // Maximumsnorm des Abstandes zweier Loesungen
  double t = 0.0;          // fuer die Zeitmessung

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

  // Initialisiere A, b und x
  init(x, b, A, n);

  it=0;

  // Beginn Zeitmessung
  struct timeval timer;
  reset_timer(&timer);

  while(!tolerance_reached(x,b,A,n) && it<((long) (maxIt))){

    jacobi(y,x,b,A,n);

    if(tolerance_reached(y,b,A,n))
      break;
    ++it;

    jacobi(x,y,b,A,n);

    ++it;
  }

  t = get_timer(timer);
  fprintf(stdout, "Iter: %d, t: %f s, Norm(res): %.16e\n", it, t, max_norm(x,b,A,n));

  // release memory
  free(x);
  free(y);
  free(b);
  free(A);

  return (0);
}
