#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BLACK 0
#define RED 1

int color(int index, int nx)
{
	floor(index/nx)+(index%nx)%2;
}

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
void init(double* x, double* b, double* A, int nx, int ny, double* hx)
{
  long i;
  long n;
  n=nx*ny;
  for (i=0; i<n; i++)
  {
    x[i] = 0.0;
    b[i] = i;
    if(i>0)A[i*n+i-1] = -1.0;
    if(i<n-1)A[i*n+i+1] = -1.0;
    if(i<n-nx-1)A[i*n+i+nx] = -1.0;
    if(i>nx)A[i*n+i-nx] = -1.0;
    A[i*n+i] = 4.0;
  }
  int c;
  for(c=0;c<2*nx;c++){
	if(c<nx)hx[c]=2.0;
	if(c>=nx)hx[(nx*ny)+nx+c]=2.0;
  }
}

// Testet ob die Abbruchbedingung erfuellt ist
int tolerance_reached(const double* x, const double* b, const double* A, const long n){
  const double eps = 1e-5; // Abbruchbedingung
  double diff;

  diff = max_norm(x,b,A,n);
  if(0)
    fprintf(stdout, "Norm(res): %.16e\n", diff);

  return diff < eps;
}

void iteration(double *u, double *b, int nx, int ny, int c, double* hx, int l, int rank)
{
	int i;
	//Hilfsvektor füllen
	for(i=nx-1;i<nx*ny+nx;i++){
		hx[i]=u[i];
	}
	for(i=rank*l+nx-1;i<(rank+1)*l+nx;i++){
		if(color(i,nx)==c){
			u[i-nx+1] = 0.25*(b[i]+hx[i-1]+hx[i+1]+hx[i-nx]+hx[i+nx]);
		}
	}
}

int main(int argc, char* argv[])
{
	int nx;
	int ny;
	int rank;
	int psize;
	const int maxIt = 10000;
	int it;
	int xinh;
	double* x;
	double* b;
	double* A;
	double* hx;
	double starttime,endtime;

	nx = atoi(argv[1]);
	ny = atoi(argv[2]);
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&psize);
	xinh=(nx*ny)/psize;
	if(rank==psize-1 && (nx*ny)%psize!=0)xinh+=(nx*ny)%psize;

	x=(double *)malloc(sizeof(double)*nx*ny);
	hx=(double *)malloc(sizeof(double)*nx*ny*2*nx);
	b=(double *)malloc(sizeof(double)*nx*ny);
	A=(double *)malloc(sizeof(double)*nx*nx*ny*ny);
	init(x,b,A,nx,ny,hx);

	if(!rank){
		starttime = MPI_Wtime();
	}

	while(!tolerance_reached(x,b,A,(nx*ny)) && it<maxIt){
		iteration(x,b,nx,ny,BLACK,hx,xinh,rank);
//		MPI_Barrier(MPI_COMM_WORLD);
		//x-sync
		int c;
		for(c=0;c<psize;c++){
			MPI_Bcast(&x[c*xinh],xinh,MPI_DOUBLE,c,MPI_COMM_WORLD);
		}
//		MPI_Barrier(MPI_COMM_WORLD);
		iteration(x,b,nx,ny,RED,hx,xinh,rank);
//		MPI_Barrier(MPI_COMM_WORLD);
		//x-sync
		for(c=0;c<psize;c++){
			MPI_Bcast(&x[c*xinh],xinh,MPI_DOUBLE,c,MPI_COMM_WORLD);
		}
//		MPI_Barrier(MPI_COMM_WORLD);
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
