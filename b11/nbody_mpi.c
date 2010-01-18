#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "io_mpi.h"
#include "generate_mpi.h"
#include "stopwatch.h"

/* avoid *[]-confusion */
typedef double double3[3];

/*
 * Global constants
 */

/* set gravity constant to one */
const double G = 1.0;
/* The epsilon^2 for the softening */
const double epsilon2 = 1E-14;

/*for tiled version*/
const int B=50;
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)>(y)) ? (x) : (y))

/**
 * Compute acceleration
 *
 * \param[in]  n Number of bodies
 * \param[in]  x Positions of the bodies
 * \param[in]  m Masses of the bodies
 * \param[out] a The computed accelerations of the bodies
 *
 * The inner loop executes 26 floating point operations and is executed
 * ((n*n)-n)/2 times, resulting in 13*n*(n-1) floating point operations total
 * per call to acceleration().
 */
void acceleration_mpi (int n, double3 x[], double m[], double3 a[])
{
  int I,J,i,j,iend,jend;
  double d0,d1,d2,r,r2,factor;
  double X[B][3], A[B][3];
  double XO[B][3], MO[B];

  /* compute acceleration *not*exploiting symmetry */
  for (I=0; I<n; I+=B) {
    /* in case n is not a multiple of B */
    iend = MIN(B,n-I);

    /* load to avoid aliasing */
    memcpy(X, x+I, sizeof(double)*3*iend);
    memcpy(A, a+I, sizeof(double)*3*iend);

    for (J=0; J<n; J+=B) {
      /* in case n is not a multiple of B */
      jend = MIN(B,n-J);

      /* load to avoid aliasing */
      memcpy(XO, x+J, sizeof(double)*3*jend);
      memcpy(MO, m+J, sizeof(double)*jend);

      /* compute locally */
      for (i=0; i<iend; i++) {
        for (j=0; j<jend; j++) {
          d0 = XO[j][0]-X[i][0]; d1 = XO[j][1]-X[i][1]; d2 = XO[j][2]-X[i][2];
          r2 = d0*d0 + d1*d1 + d2*d2 + epsilon2;
          r = sqrt(r2);
          factor = MO[j]*G/(r*r2);
          A[i][0] += factor*d0; A[i][1] += factor*d1; A[i][2] += factor*d2;
        }
      }
    }
    /* store result */
    memcpy(a+I, A, sizeof(double)*3*iend);
  }
}

/**
 * Do one time step
 *
 * \param[in]     n    Number of bodies
 * \param[in]     dt   Time step size
 * \param[in,out] x    IN: old positions, OUT: new positions
 * \param[in,out] v    IN: old velocities, OUT: new velocities
 * \param[in]     m    Masses
 * \param[in,out] a    IN: old accelerations, OUT: new accelerations
 * \param         aold Buffer to hold the old accelerations.  Input and output
 *                     contents is not meaningful.
 *
 * The "update position" and "update velocity" execute 12 floating point
 * operations per iteration each, resulting in 24*n floating point ops
 * together.  Add to it the 3 flops to precompute time step related values and
 * the 13*n*(n-1) flops from acceleration() to arrive at the total of
 * 13*n*(n-1)+24*n+2 floating point operations per call to leapfrog().
 */
void leapfrog (int n, double dt, double3 x[], double3 v[], double m[],
               double3 a[], double3 aold[])
{
  int i;
  /* precompute some values from the time step */
  double dt2 = dt*dt*0.5;
  double dthalf = dt*0.5;

  /* update position */
  for (i=0; i<n; i++) {
    x[i][0] += dt*v[i][0] + dt2*a[i][0];
    x[i][1] += dt*v[i][1] + dt2*a[i][1];
    x[i][2] += dt*v[i][2] + dt2*a[i][2];
  }

  /* save and clear acceleration */
  for (i=0; i<n; i++) {
    aold[i][0] = a[i][0];
    aold[i][1] = a[i][1];
    aold[i][2] = a[i][2];
    a[i][0] = a[i][1] = a[i][2] = 0.0;
  }

  /* compute new acceleration */
  acceleration_mpi(n,x,m,a);

  /* update velocity */
  for (i=0; i<n; i++) {
    v[i][0] += dthalf*aold[i][0] + dthalf*a[i][0];
    v[i][1] += dthalf*aold[i][1] + dthalf*a[i][1];
    v[i][2] += dthalf*aold[i][2] + dthalf*a[i][2];
  }
}

/**
 * main program
 */
int main (int argc, char** argv)
{
  /* hold positions, velocities, masses, and accelerations for the bodies */
  double3 *x;
  double3 *v;
  double *m;
  double3 *a;
  /* temporary buffer given to leapfrog() so it can store the old acceleration
     values there */
  double3 *temp;
  /* loop counter */
  int i;
  /* count the time steps done */
  int k = 0;
  /* base name for the output files */
  char base[]="nbody_mpi";
  /* buffer to construct the whole output filename in */
  char name[256];
  /* handle for the output file */
  FILE *file;
  /* accumulated simulation time and time step size.  These will be modified
     when reading initial values from a file */
  double t = 0, dt = 1E-3;
  /* for measuring the elapsed time and calculating the MFLOP rate */
  double start,stop,elapsed,flop;
  /* N - number of bodies (read from commandline), n - per processor */
  int n, N;
  /* number of time steps to do (read from commandline) */
  int timesteps;
  /* when to measure MFLOP rate and write output files: every mod time steps
     (read from commandline) */
  int mod;

  int rank, procs;

  // MPI initialisieren
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  /**
   * parse commandline arguments
   */
  if (argc!=4) {
    printf("usage: nbody_mpi <nbodies> <timesteps> <every>\n");
    return 1;
  }
  sscanf(argv[1],"%d",&N);
  sscanf(argv[2],"%d",&timesteps);
  sscanf(argv[3],"%d",&mod);

  if(N%procs != 0) {
    fprintf(stderr,
            "<nbodies> must be a multiple of the number of processes %d\n",
            procs);
    return 1;
  }
  n = N/procs;

  x = calloc(n,sizeof(double3));
  v = calloc(n,sizeof(double3));
  m = calloc(n,sizeof(double));
  plummer_mpi(n,17,x,v,m);

  a = calloc(n,sizeof(double3));
  temp = calloc(n,sizeof(double3));
  for (i=0; i<n; i++)
    a[i][0] = a[i][1] = a[i][2] = 0.0;
  acceleration_mpi(n,x,m,a);

  k = 0;
  t = 0.0;
  dt = 1E-3;
  if(rank == 0) {
    sprintf(name,"%s_%06d.vtk",base,k);
    printf("writing %s \n",name);
    file = fopen(name,"w");
  }
  write_vtk_file_double_mpi(file,n,x,v,m,t,dt);
  if(rank == 0)
    fclose(file);
  start = get_time();

  for (k=1; k<=timesteps; k++) {
    leapfrog(n,dt,x,v,m,a,temp);
    t += dt;
    if (k%mod==0) {
      stop = get_time();
      elapsed = stop-start;
      flop = mod*(19.0*N*N+24.0*N);
      if(rank == 0) {
        printf("%g seconds for %g ops = %g MFLOPS\n",elapsed,flop,flop/elapsed/1E6);

        sprintf(name,"%s_%06d.vtk",base,k/mod);
        printf("writing %s \n",name);
        file = fopen(name,"w");
      }
      write_vtk_file_double_mpi(file,n,x,v,m,t,dt);
      if(rank == 0)
        fclose(file);

      start = get_time();
    }
  }

  /* Don't bother freeing the allocated buffers since the program is going to
     end and the operating system will reclaim all of its memory anyway */

  // MPI geordnet beenden
  MPI_Finalize();

  return 0;
}
