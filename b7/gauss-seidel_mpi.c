#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

// Testet ob die Abbruchbedingung erfuellt ist
bool tolerance_reached(const double* x, const double* b, const double* A, const long n){
  const double eps = 1e-5; // Abbruchbedingung
  double diff;

  diff = max_norm(x,b,A,n);
  if(false)
    fprintf(stdout, "Norm(res): %.16e\n", diff);

  return diff < eps;
}

void iteration(double *u, double *x, double *b, long n)
{
	int i;
	for(i=0;i<n*n;i++){
		if(i<=nr){
			u[i] = 0.25*(b[i]+u[i-1]+u[i+1]+u[i-nx]+u[i+nx]);
		} else {
			u[i] = 0.25*(b[i]+u[i-1]+u[i+1]+u[i-nx]+u[i+nx]);
		}
	}
}

int main(int argc, char* argv[])
{
	long n;
	int rank;
	int psize;
	const int maxIt = 10000;
	int it;
	int xinh;
	double* x;
	double* b;
	double* A;
	double starttime,endtime;

	if(argc!=2){
		printf("Usage: %s <x>\n",argv[0]);
		exit(1);
	}
	n = atol(argv[1]);
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&psize);

	if(!rank){
		x=(double *)malloc(sizeof(double)*n);
		b=(double *)malloc(sizeof(double)*n);
		A=(double *)malloc(sizeof(double)*n*n);
		init(x,b,A,n);
		starttime = MPI_Wtime();
	}

	//groesse von x bestimmen
	xinh = n/psize;
	if(x%n!=0 && rank==psize-1)xinh+=x%n;

	while(!tolerance_reached(x,b,A,n) && it<maxIt){
		//x-sync
		for(c=0;c<psize;c++){
			MPI_Bcast(&x,xinh,MPI_DOUBLE,rank,MPI_COMM_WORLD);
		}
		iteration(&A,&x,&b,n);
		++it;
	}

	if(!rank){
		endtime = MPI_Wtime();
		printf("%f\n",endtime-starttime);
	}



	MPI_Finalize();
	free(x);
	free(b);
	free(A);
}
