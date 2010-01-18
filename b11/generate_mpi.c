#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double double3[3];

void cube_mpi (int n, long int seed, double size, double m0, double mdelta,  
               double x[][3], double v[][3], double m[])
{
  int i, j;
  double s[3];
  double t[3];
  double s_all[3];
  double t_all[3];
  double M;
  double M_all;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (seed!=0) srand48(seed);

  // loop to skip mass points of other processes
  for (j=0; j<rank+1; j++)
    {
      s[0] = 0.0; s[1] = 0.0; s[2] = 0.0;
      t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;
      M = 0.0;
      for (i=0; i<n; i++)
        {
          x[i][0] = drand48()*size;
          x[i][1] = drand48()*size;
          x[i][2] = drand48()*size;
          v[i][0] = 0.0;
          v[i][1] = 0.0;
          v[i][2] = 0.0;
          m[i] = m0 + (drand48()-0.5)*2.0*mdelta;
          s[0] += m[i]*x[i][0];
          s[1] += m[i]*x[i][1];
          s[2] += m[i]*x[i][2];
          t[0] += m[i]*v[i][0];
          t[1] += m[i]*v[i][1];
          t[2] += m[i]*v[i][2];
          M += m[i];
        }
    }
  MPI_Allreduce(s, s_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(t, t_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&M, &M_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  s_all[0] /= M_all; s_all[1] /= M_all; s_all[2] /= M_all;
  t_all[0] /= M_all; t_all[1] /= M_all; t_all[2] /= M_all;
  if (rank == 0)
    printf("center of mass: %g %g %g\n",s_all[0],s_all[1],s_all[2]);
  for (i=0; i<n; i++)
    {
      x[i][0] -= s_all[0]; 
      x[i][1] -= s_all[1]; 
      x[i][2] -= s_all[2]; 

      v[i][0] -= t_all[0]; 
      v[i][1] -= t_all[1]; 
      v[i][2] -= t_all[2]; 
    }
}

void plummer_mpi (int n, long int seed,
                  double x[][3], double v[][3], double m[])
{
  /* This is a copy from 
     Pit Hut, Jun Makino: The Art of Computational Science, 
     The Kali Code, vol. 5. Initial Conditions: Plummer's Model.
  */ 
  int i, j;
  double radius,theta,phi,X,Y,velocity,maxr;
  const double Pi = 3.141592653589793238462643383279;
  double s[3];
  double t[3];
  double s_all[3];
  double t_all[3];
  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if (seed!=0) srand48(seed);

  // loop to skip mass points of other processes
  for (j=0; j<rank+1; j++)
    {
      s[0] = 0.0; s[1] = 0.0; s[2] = 0.0;
      t[0] = 0.0; t[1] = 0.0; t[2] = 0.0;
      maxr = -1.0;
      for (i=0; i<n; i++)
        {
          m[i] = 1.0/(n*procs);
          radius = 1.0/sqrt(pow(drand48(),-2.0/3.0)-1.0);
          if (radius>maxr) maxr=radius;
          theta = acos(-1.0+drand48()*2.0);
          phi = drand48()*2.0*Pi;
          x[i][0] = radius*sin(theta)*cos(phi);
          x[i][1] = radius*sin(theta)*sin(phi);
          x[i][2] = radius*cos(theta);
          s[0] += m[i]*x[i][0];
          s[1] += m[i]*x[i][1];
          s[2] += m[i]*x[i][2];
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
          v[i][0] = velocity*sin(theta)*cos(phi);
          v[i][1] = velocity*sin(theta)*sin(phi);
          v[i][2] = velocity*cos(theta);
          t[0] += m[i]*v[i][0];
          t[1] += m[i]*v[i][1];
          t[2] += m[i]*v[i][2];
        }
    }
  MPI_Allreduce(s, s_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(t, t_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf("center of mass: %g %g %g\n",s_all[0],s_all[1],s_all[2]);
  for (i=0; i<n; i++)
    {
      x[i][0] -= s_all[0]; 
      x[i][1] -= s_all[1]; 
      x[i][2] -= s_all[2]; 

      v[i][0] -= t_all[0]; 
      v[i][1] -= t_all[1]; 
      v[i][2] -= t_all[2]; 
    }
  s[0]=s[1]=s[2]=0.0;
  for (i=0; i<n; i++)
    {
      s[0] += m[i]*x[i][0];
      s[1] += m[i]*x[i][1];
      s[2] += m[i]*x[i][2];
    }
  MPI_Allreduce(s, s_all, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf("new center of mass: %g %g %g\n",s_all[0],s_all[1],s_all[2]);
}
