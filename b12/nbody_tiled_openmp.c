#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "io_vanilla.h"
#include "generate_vanilla.h"
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
/* Make the tile size a constant to allow better loop-unrolling by the
   compiler */
const int B = 1;

#define acceleration acceleration_tiled

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
void acceleration_tiled (int n, double3 x[], double m[], double3 a[])
{
  int i,j,I,J,iend,jend,jstart;
  double d0,d1,d2,r,r2,factori,factorj;
  double invfact;

  /* compute acceleration exploiting symmetry */
  for (I=0; I<n; I+=B) {
    iend = MIN(I+B, n);
    for (J=I+1; J<n; J+=B) {
      jend = MIN(J+B,n);

      for (i=I; i<iend; i++) {
        jstart = MAX(i+1, J);

        for (j=jstart; j<jend; j++) {
          /* compute distance vector */
          d0 = x[j][0]-x[i][0];
          d1 = x[j][1]-x[i][1];
          d2 = x[j][2]-x[i][2];

          /* compute (scalar valued) distance and distance^2, with softening */
          r2 = d0*d0 + d1*d1 + d2*d2 + epsilon2;
          r = sqrt(r2);

          /* precompute the factors in the acceleration equation */
          invfact = G/(r*r2);
          factori = m[i]*invfact;
          factorj = m[j]*invfact;

          /* compute acceleration of body i due to body j */
          a[i][0] += factorj*d0;
          a[i][1] += factorj*d1;
          a[i][2] += factorj*d2;

          /* compute acceleration of body j due to body i */
          a[j][0] -= factori*d0;
          a[j][1] -= factori*d1;
          a[j][2] -= factori*d2;
        }
      }
    }
  }
}

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
void acceleration_tiled_local (int n, double3 x[], double m[], double3 a[])
{
  int i,j,I,J,iend,jend,jstart;
  double d0,d1,d2,r,r2,factori,factorj;
  double invfact;
  double3 xli[B], xlj[B];
  double mli[B], mlj[B];
  double3 ali[B], alj[B];

  /* compute acceleration exploiting symmetry */
  for (I=0; I<n; I+=B) {
    iend = MIN(B, n-I);
    /* copy into local variables to avoid aliasing */
    for(i = 0; i < iend; ++i) {
      xli[i][0] = x[I+i][0]; xli[i][1] = x[I+i][1]; xli[i][2] = x[I+i][2];
      mli[i] = m[I+i];
      ali[i][0] = ali[i][1] = ali[i][2] = 0;
    }
    for (J=I+1; J<n; J+=B) {
      jend = MIN(B,n-J);
      /* copy into local variables to avoid aliasing */
      for(j = 0; j < jend; ++j) {
        xlj[j][0] = x[J+j][0]; xlj[j][1] = x[J+j][1]; xlj[j][2] = x[J+j][2];
        mlj[j] = m[J+j];
        alj[j][0] = alj[j][1] = alj[j][2] = 0;
      }

      for (i=0; i<iend; i++) {
        jstart = MAX(J,I+i+1)-J;

        for (j=jstart; j<jend; j++) {
          /* compute distance vector */
          d0 = xlj[j][0]-xli[i][0];
          d1 = xlj[j][1]-xli[i][1];
          d2 = xlj[j][2]-xli[i][2];

          /* compute (scalar valued) distance and distance^2, with softening */
          r2 = d0*d0 + d1*d1 + d2*d2 + epsilon2;
          r = sqrt(r2);

          /* precompute the factors in the acceleration equation */
          invfact = G/(r*r2);
          factori = mli[i]*invfact;
          factorj = mlj[j]*invfact;

          /* compute acceleration of body i due to body j */
          ali[i][0] += factorj*d0;
          ali[i][1] += factorj*d1;
          ali[i][2] += factorj*d2;

          /* compute acceleration of body j due to body i */
          alj[j][0] -= factori*d0;
          alj[j][1] -= factori*d1;
          alj[j][2] -= factori*d2;
        }
      }

      /* accumulate back to global variables */
      for(j = 0; j < jend; ++j) {
        a[J+j][0] += alj[j][0];
        a[J+j][1] += alj[j][1];
        a[J+j][2] += alj[j][2];
      }
    }
    for(i = 0; i < iend; ++i) {
      a[I+i][0] += ali[i][0];
      a[I+i][1] += ali[i][1];
      a[I+i][2] += ali[i][2];
    }
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
  for (i=0; i<n; i++)
    {
      x[i][0] += dt*v[i][0] + dt2*a[i][0];
      x[i][1] += dt*v[i][1] + dt2*a[i][1];
      x[i][2] += dt*v[i][2] + dt2*a[i][2];
    }

  /* save and clear acceleration */
  for (i=0; i<n; i++)
    {
      aold[i][0] = a[i][0];
      aold[i][1] = a[i][1];
      aold[i][2] = a[i][2];
      a[i][0] = a[i][1] = a[i][2] = 0.0;
    }

  /* compute new acceleration */
  acceleration(n,x,m,a);

  /* update velocity */
  for (i=0; i<n; i++)
    {
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
  char base[]="test1";
  /* buffer to construct the whole output filename in */
  char name[256];
  /* handle for the output file */
  FILE *file;
  /* accumulated simulation time and time step size.  These will be modified
     when reading initial values from a file */
  double t = 0, dt = 1E-3;
  /* for measuring the elapsed time and calculating the MFLOP rate */
  double start,stop,elapsed,flop;
  /* number of bodies (read from commandline) */
  int n;
  /* number of time steps to do (read from commandline) */
  int timesteps;
  /* when to measure MFLOP rate and write output files: every mod time steps
     (read from commandline) */
  int mod;
  /* number of OpenMP threads to use (read from commandline) */
  int numthreads;

  /**
   * parse commandline arguments
   */
  if (argc!=5) {
    printf("usage: nbody_vanilla <nbodies> <timesteps> <every> <numthreads>\n");
    return 1;
  }
  sscanf(argv[1],"%d",&n);
  sscanf(argv[2],"%d",&timesteps);
  sscanf(argv[3],"%d",&mod);
  sscanf(argv[4],"%d",&numthreads);

#ifdef _OPENMP
  printf("number of processors: %d\n",omp_get_num_procs());
  omp_set_num_threads(numthreads);
  printf("number of threads used: %d\n",omp_get_max_threads());
#else
  printf("not compiled with OpenMP!\n");
  numthreads = 1;
#endif

#if 0
  /* get initial value from a data file */
  file = fopen("plummer10000.vtk","r");
  if (file==NULL) {
    printf("could not open file --- aborting\n");
    return 1;
  }
  /* peek into the file to find out the number of bodies */
  n = get_vtk_numbodies(file);
  rewind(file);
  /* allocate memory accordingly */
  x = calloc(n,sizeof(double3));
  v = calloc(n,sizeof(double3));
  m = calloc(n,sizeof(double));
  /* read everything */
  read_vtk_file_double(file,n,x,v,m,&t,&dt);
  fclose(file);
  printf("loaded %d bodies\n",n);
#else
  /* get initial values from one of the generator functions */
  x = calloc(n,sizeof(double3));
  v = calloc(n,sizeof(double3));
  m = calloc(n,sizeof(double));
  plummer(n,17,x,v,m);
#endif

  /**
   * allocate and initialise the acceleration and the buffer for leapfrog
   */
  a = calloc(n,sizeof(double3));
  for (i=0; i<n; i++)
    a[i][0] = a[i][1] = a[i][2] = 0.0;
  acceleration(n,x,m,a);
  temp = calloc(n,sizeof(double3));

  /* write out the initial values */
  sprintf(name,"%s_%06d.vtk",base,k/mod);
  printf("writing %s \n",name);
  file = fopen(name,"w");
  write_vtk_file_double(file,n,x,v,m,t,dt);
  fclose(file);

  /* start the time measurement */
  start = get_time();
  for (k=1; k<timesteps; k++) {
    /* do one timestep and increment elapsed simulated time */
    leapfrog(n,dt,x,v,m,a,temp);
    t += dt;
    if (k%mod==0) {
      /* write statistics */
      stop = get_time();
      elapsed = stop-start;
      /* 13*n*(n-1)+24*n+3 FLOP from leaprog(), 1 from the t+=dt above */
      flop = mod*(13.0*n*(n-1.0)+24.0*n+4.0);
      printf("%g seconds for %g ops = %g MFLOPS\n",
             elapsed,flop,flop/elapsed/1E6);

      /* write output file */
      sprintf(name,"%s_%06d.vtk",base,k/mod);
      printf("writing %s \n",name);
      file = fopen(name,"w");
      write_vtk_file_double(file,n,x,v,m,t,dt);
      fclose(file);

      /* restart the timer */
      start = get_time();
    }
  }

  /* Don't bother freeing the allocated buffers since the program is going to
     end and the operating system will reclaim all of its memory anyway */

  return 0;
}
