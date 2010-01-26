/* Necessary changes to parallelize the vanilla version with cuda

   (1) copy all subroutines for generating, I/O and time measurement to here
       (mainly necessary to avoid complicating the Makefile).
   (2) make it compile with C++ compiler (as nvcc uses that)
   (3) replace double3 by float4, replace all x[][3] by float4*
   (4) store mass in fourth component of position
   (5) use nvcc
   (6) float4 is struct in CUDA
   (7) remove symmetry optimization

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>

// typedef float float4[4]; // provided by NVIDIA

/*const double gamma = 6.674E-11;*/
const float G = 1.0;
const float epsilon2 = 1E-5;

// number of threads in one block
#define BLOCKSIZE 128

// kernel to be executed on the GPU
__global__ void acceleration_kernel (int jstart, float4 *x, float4 *a)
{
	extern __shared__ float4 xcol[];
    int tid = threadIdx.x;
	int i = blockIdx.x*blockDim.x+tid;

    // load column positions an mass in shared memory, each thread loads one element 
	xcol[tid] = x[jstart+tid];
	__syncthreads();

	// compute contribution to body i
	float4 xi = x[i];  // load own position
	float4 ai = a[i];  // place to accumulate result for own body
  	float3 d;
	float r,r2,factorj,invfact;
	for (int j=0; j<blockDim.x; j++)
      {
		d.x = xcol[j].x-xi.x;
		d.y = xcol[j].y-xi.y;
		d.z = xcol[j].z-xi.z;
		r2 = d.x*d.x + d.y*d.y + d.z*d.z + epsilon2;
		r = sqrtf(r2);
		invfact = G/(r*r2);
		factorj = xcol[j].w*invfact;
		ai.x += factorj*d.x;
		ai.y += factorj*d.y;
		ai.z += factorj*d.z;
      }

    // write back result
    a[i] = ai;
}


void acceleration (int n, float4 *x, float4 *a)
{
    int size = n*sizeof(float4);

    // allocate x in global memory on the device
    float4 *xd;
    cudaMalloc( (void**) &xd, size ); // allocate memory on device
    cudaMemcpy(xd,x,size,cudaMemcpyHostToDevice); // copy x to device
    if( cudaGetLastError() != cudaSuccess) 
    {
        fprintf(stderr,"error in memcpy\n");
        exit(-1);
    }                         

    // allocate a in global memory on the device
	// note: a has already been cleared
    float4 *ad;
    cudaMalloc( (void**) &ad, size ); // allocate memory on device
    cudaMemcpy(ad,a,size,cudaMemcpyHostToDevice); // copy a to device
    if( cudaGetLastError() != cudaSuccess) 
    {
        fprintf(stderr,"error in memcpy\n");
        exit(-1);
    }                         

	// determine block and grid size
	dim3 dimBlock(BLOCKSIZE);   // use BLOCKSIZE threads in one block
    dim3 dimGrid(n/BLOCKSIZE);  // use enough blocks to reach n; n is a multiple of BLOC

	// process chunks of rows
	for (int j=0; j<n; j+=BLOCKSIZE)
    {
        acceleration_kernel<<<dimGrid,dimBlock,BLOCKSIZE*sizeof(float4)>>>(j,xd,ad);
  		cudaThreadSynchronize();
    }

	// read result
	cudaMemcpy(a,ad,size,cudaMemcpyDeviceToHost);
    if( cudaGetLastError() != cudaSuccess) 
    {
        fprintf(stderr,"error in memcpy\n");
        exit(-1);
    }                         

	// free memory on device
	cudaFree(xd);
	cudaFree(ad);
}

void acceleration_cpu (int n, float4 *x, float4 *a)
{
  int i,j;
  float d0,d1,d2,r,r2,factorj,invfact;

  /* compute acceleration exploiting symmetry */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      {
		d0 = x[j].x-x[i].x;
		d1 = x[j].y-x[i].y;
		d2 = x[j].z-x[i].z;
		r2 = d0*d0 + d1*d1 + d2*d2 + epsilon2;
		r = sqrt(r2);
		invfact = G/(r*r2);
		factorj = x[j].w*invfact;
		a[i].x += factorj*d0;
		a[i].y += factorj*d1;
		a[i].z += factorj*d2;
      }
}

void leapfrog (int n, float dt, float4 *x, float4 *v, float4 *a, float4 *aold)
{
  int i;
  float dt2 = dt*dt*0.5;
  float dthalf = dt*0.5;

  /* update position */
  for (i=0; i<n; i++)
    {
      x[i].x += dt*v[i].x + dt2*a[i].x;
      x[i].y += dt*v[i].y + dt2*a[i].y;
      x[i].z += dt*v[i].z + dt2*a[i].z;
    }

  /* save and clear acceleration */
  for (i=0; i<n; i++)
    {
      aold[i].x = a[i].x;
      aold[i].y = a[i].y;
      aold[i].z = a[i].z;
      a[i].x = a[i].y = a[i].z = 0.0;
    }
  
  /* compute new acceleration */
  acceleration(n,x,a);

  /* update velocity */
  for (i=0; i<n; i++)
    {
      v[i].x += dthalf*aold[i].x + dthalf*a[i].x;
      v[i].y += dthalf*aold[i].y + dthalf*a[i].y;
      v[i].z += dthalf*aold[i].z + dthalf*a[i].z;
    }
}


void write_vtk_file_float (FILE *f, int n, float4 *x, float4 *v,
			    float t, float dt)
{
  int i;

  /* header */
  fprintf(f,"%s\n","# vtk DataFile Version 1.0");
  fprintf(f,"NBODY %22.16g %22.16g\n",t,dt);
  fprintf(f,"%s\n","ASCII");

  /* points */
  fprintf(f,"%s\n","DATASET POLYDATA");
  fprintf(f,"%s %d %s\n","POINTS",n,"float");
  for (i=0; i<n; i++)
    fprintf(f,"%22.16g %22.16g %22.16g\n",x[i].x,x[i].y,x[i].z);

  /* vertices */
  fprintf(f,"%s %d %d\n","VERTICES",n,2*n);
  for (i=0; i<n; i++)
    fprintf(f,"%d %d\n",1,i);

  /* scalar data fields*/
  fprintf(f,"%s %d\n","POINT_DATA",n);
  fprintf(f,"%s\n","SCALARS mass float");
  fprintf(f,"%s\n","LOOKUP_TABLE default");
  for (i=0; i<n; i++)
    fprintf(f,"%22.16g\n",x[i].w);

  /* vector data */
  fprintf(f,"%s\n","VECTORS velocity float");
  for (i=0; i<n; i++)
    fprintf(f,"%22.16g %22.16g %22.16g\n",v[i].x,v[i].y,v[i].z);
}
    
void cube (int n, long int seed, float size, float m0, float mdelta,  
		float4 *x, float4 *v)
{
  int i;
  float3 s;
  float3 t;
  float M = 0.0;

  s.x = s.y = s.z = 0.0;
  t.x = t.y = t.z = 0.0;

  if (seed!=0) srand48(seed);
  for (i=0; i<n; i++)
    {
      x[i].x = drand48()*size;
      x[i].y = drand48()*size;
      x[i].z = drand48()*size;
      v[i].x = 0.0;
      v[i].y = 0.0;
      v[i].z = 0.0;
      x[i].w = m0 + (drand48()-0.5)*2.0*mdelta;
      s.x += x[i].w*x[i].x;
      s.y += x[i].w*x[i].y;
      s.z += x[i].w*x[i].z;
      t.x += x[i].w*v[i].x;
      t.y += x[i].w*v[i].y;
      t.z += x[i].w*v[i].z;
      M += x[i].w;
    }
  printf("center of mass: %g %g %g\n",s.x,s.y,s.z);
  for (i=0; i<n; i++)
    {
      x[i].x -= s.x; 
      x[i].y -= s.y; 
      x[i].z -= s.z; 

      v[i].x -= t.x; 
      v[i].y -= t.y; 
      v[i].z -= t.z; 
    }
}

void plummer (int n, long int seed,
			  float4 *x, float4 *v)
{
  /* This is a copy from 
     Pit Hut, Jun Makino: The Art of Computational Science, 
     The Kali Code, vol. 5. Initial Conditions: Plummer's Model.
  */ 
  int i;
  float radius,theta,phi,X,Y,velocity,maxr=-1.0;
  const float Pi = 3.141592653589793238462643383279;
  float3 s;
  float3 t;
  s.x = s.y = s.z = 0.0;
  t.x = t.y = t.z = 0.0;

  if (seed!=0) srand48(seed);
  for (i=0; i<n; i++)
    {
      x[i].w = 1.0/n;
      radius = 1.0/sqrt(pow(drand48(),-2.0/3.0)-1.0);
      if (radius>maxr) maxr=radius;
      theta = acos(-1.0+drand48()*2.0);
      phi = drand48()*2.0*Pi;
      x[i].x = radius*sin(theta)*cos(phi);
      x[i].y = radius*sin(theta)*sin(phi);
      x[i].z = radius*cos(theta);
      s.x += x[i].w*x[i].x;
      s.y += x[i].w*x[i].y;
      s.z += x[i].w*x[i].z;
      X = 0.0;
      Y = 0.1;
      while (Y>X*X*pow(1.0-X*X,3.5))
	{
	  X = drand48();
	  Y = drand48()*0.1;
	}
      velocity = X*sqrt(2.0)*pow(1.0+radius*radius,-0.25);
      theta = acos(-1.0+drand48()*2.0);
      phi = drand48()*2.0*Pi;
      v[i].x = velocity*sin(theta)*cos(phi);
      v[i].y = velocity*sin(theta)*sin(phi);
      v[i].z = velocity*cos(theta);
      t.x += x[i].w*v[i].x;
      t.y += x[i].w*v[i].y;
      t.z += x[i].w*v[i].z;
    }
  printf("center of mass: %g %g %g\n",s.x,s.y,s.z);
  for (i=0; i<n; i++)
    {
      x[i].x -= s.x; 
      x[i].y -= s.y; 
      x[i].z -= s.z; 

      v[i].x -= t.x; 
      v[i].y -= t.y; 
      v[i].z -= t.z; 
    }
  s.x=s.y=s.z=0.0;
  for (i=0; i<n; i++)
    {
      s.x += x[i].w*x[i].x;
      s.y += x[i].w*x[i].y;
      s.z += x[i].w*x[i].z;
    }
  printf("new center of mass: %g %g %g\n",s.x,s.y,s.z);
  printf("maximum radius: %g\n",maxr);
}

float get_time ()
{
/*   struct rusage ru; */
/*   struct timeval cstop; */
  
/*   getrusage(RUSAGE_SELF, &ru); */
/*   cstop = ru.ru_utime; */
/*   return 1.0*cstop.tv_sec + cstop.tv_usec/1000000.0; */
/*   return clock()/CLOCKS_PER_SEC; */
  return omp_get_wtime();
}

int main (int argc, char** argv)
{
  int n;
  float4 *x;
  float4 *v;
  float4 *a;
  float4 *temp;
  int i,k,mod;
  char base[]="test1";
  char name[256];
  FILE *file;
  float t,dt;
  float start,stop;
  float elapsed,flop;
  int timesteps;

  if (argc!=4)
	{
	  printf("usage: nbody_cuda <nbodies> <timesteps> <every>\n");
	  return 1;
	}
  sscanf(argv[1],"%d",&n);
  sscanf(argv[2],"%d",&timesteps);
  sscanf(argv[3],"%d",&mod);

    // determine block size and grid
    if (n%BLOCKSIZE != 0)
    {
        printf("n must be a multiple of BLOCKSIZE\n");
		exit(-1);
    }

  /* allocate memory, read data file */
  /*   file = fopen("plummer10000.vtk","r"); */
  /*   if (file==NULL) */
  /*     { */
  /*       printf("could not open file --- aborting\n"); */
  /*       return 1; */
  /*     } */
  /*   n = get_vtk_numbodies(file); */
  /*   rewind(file); */
  /*   x = calloc(n,sizeof(float4)); */
  /*   v = calloc(n,sizeof(float4)); */
  /*   m = calloc(n,sizeof(float)); */
  /*   read_vtk_file_float(file,n,x,v,m,&t,&dt); */
  /*   fclose(file); */
  /*   printf("loaded %d bodies\n",n); */
  x = (float4*) calloc(n,sizeof(float4));
  v = (float4*)calloc(n,sizeof(float4));
  plummer(n,17,x,v);

  a = (float4*) calloc(n,sizeof(float4));
  temp = (float4*) calloc(n,sizeof(float4));
  for (i=0; i<n; i++)
    a[i].x = a[i].y = a[i].z = 0.0;
  acceleration(n,x,a);

  k = 0;
  t = 0.0;
  dt = 1E-3;
  printf("writing %s_%06d.vtk \n",base,k);
  sprintf(name,"%s_%06d.vtk",base,k);
  file = fopen(name,"w");
  write_vtk_file_float(file,n,x,v,t,dt);
  fclose(file);
  start = get_time();

  for (k=1; k<timesteps; k++)
    {
      leapfrog(n,dt,x,v,a,temp);
      t += dt;
      if (k%mod==0)
		{
		  stop = get_time();
		  elapsed = stop-start;
		  flop = mod*(19.0*n*n+24.0*n);
		  printf("%g seconds for %g ops = %g MFLOPS\n",elapsed,flop,flop/elapsed/1E6);
		  printf("writing %s_%06d.vtk \n",base,k/mod);
		  
		  sprintf(name,"%s_%06d.vtk",base,k/mod);
		  file = fopen(name,"w");
		  write_vtk_file_float(file,n,x,v,t,dt);
		  fclose(file);
		  
		  start = get_time();
		}
    }

  return 0;
}
