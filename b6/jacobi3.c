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
void jacobi(double* xneu, const double* const xalt, const double* const b, const double* const A, long start, long width, long n)
{
  int i=0, j=0;
/*printf("xneu(0) %f\n",xneu[0]);
printf("xneu(1) %f\n",xneu[1]);
printf("xalt(0) %f\n",xalt[0]);
printf("xalt(1) %f\n",xalt[1]);
printf("b(0) %f\n",b[0]);
printf("b(1) %f\n",b[1]);
printf("A(0) %f\n",A[0]);
printf("A(0) %f\n",A[1]);
printf("A(0) %f\n",A[2]);
printf("A(0) %f\n",A[3]);*/

  for (i=start; i<start+width; ++i)
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
  xx = (double *) malloc(sizeof(double)*n);
  yy = (double *) malloc(sizeof(double)*n);

  // Beginn Zeitmessung
  struct timeval timer;
  if(rank==0)reset_timer(&timer);

  // Jacobi-Schritte
  for (it=1; it<((long) (maxIt/2)+1); it++)
  {
    // Trick um das Umkopieren x^{m} = y zu sparen:
    // 2 Jacobi-Schritte, Loesung x^(m) ist dann
    // abwechselnd in x oder y
    MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(yy,x,b,A,rank*sendcounts[0],sendcounts[rank],n);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(yy,sendcounts[rank],MPI_DOUBLE,y,sendcounts,senddispls,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//printf("rank(%d) y(0): %f\n",rank,y[0]);
	//printf("rank(%d) y(1): %f\n",rank,y[1]);
	//printf("rank(%d) y(2): %f\n",rank,y[2]);
	//printf("rank(%d) y(3): %f\n",rank,y[3]);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(y,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    jacobi(xx,y,b,A,rank*sendcounts[0],sendcounts[rank],n);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(xx,sendcounts[rank],MPI_DOUBLE,x,sendcounts,senddispls,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	//printf("rank(%d) x(0): %f\n",rank,x[0]);
        //printf("rank(%d) x(1): %f\n",rank,x[1]);
        //printf("rank(%d) x(2): %f\n",rank,x[2]);
        //printf("rank(%d) x(3): %f\n",rank,x[3]);
	
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
  free(xx);
  free(yy);
  free(b);
  free(A);
  MPI_Finalize();
  return (0);
}
