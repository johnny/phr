#include <stdio.h>

#include <mpi.h>

#include "io_mpi.h"

void write_vtk_file_double_mpi (FILE *f, int n, double x[][3], double v[][3],
                                double m[], double t, double dt)
{
  int i;
  int rank, procs;
  int p;
  double vec[n][3], sca[n];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  /* header */
  if(rank == 0) {
    fprintf(f,"%s\n","# vtk DataFile Version 1.0");
    fprintf(f,"NBODY %22.16g %22.16g\n",t,dt);
    fprintf(f,"%s\n","ASCII");
  }

  /* points */
  if(rank == 0) {
    fprintf(f,"%s\n","DATASET POLYDATA");
    fprintf(f,"%s %d %s\n","POINTS",n*procs,"float");
    for (i=0; i<n; i++)
      fprintf(f,"%22.16g %22.16g %22.16g\n",x[i][0],x[i][1],x[i][2]);
    for (p = 1; p < procs; ++p) {
      MPI_Recv(vec, n*3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (i=0; i<n; i++)
        fprintf(f,"%22.16g %22.16g %22.16g\n",vec[i][0],vec[i][1],vec[i][2]);
    }
  } else {
    MPI_Ssend(x, n*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  /* vertices */
  if(rank == 0) {
    fprintf(f,"%s %d %d\n","VERTICES",n*procs,2*n*procs);
    for (i=0; i<n*procs; i++)
      fprintf(f,"%d %d\n",1,i);
  }

  /* scalar data fields*/
  if(rank == 0) {
    fprintf(f,"%s %d\n","POINT_DATA",n*procs);
    fprintf(f,"%s\n","SCALARS mass float");
    fprintf(f,"%s\n","LOOKUP_TABLE default");
    for (i=0; i<n; i++)
      fprintf(f,"%22.16g\n",m[i]);
    for (p = 1; p < procs; ++p) {
      MPI_Recv(sca, n, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (i=0; i<n; i++)
        fprintf(f,"%22.16g\n",sca[i]);
    }
  } else {
    MPI_Ssend(m, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  /* vector data */
  if(rank == 0) {
    fprintf(f,"%s\n","VECTORS velocity float");
    for (i=0; i<n; i++)
      fprintf(f,"%22.16g %22.16g %22.16g\n",v[i][0],v[i][1],v[i][2]);
    for (p = 1; p < procs; ++p) {
      MPI_Recv(vec, n*3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (i=0; i<n; i++)
        fprintf(f,"%22.16g %22.16g %22.16g\n",vec[i][0],vec[i][1],vec[i][2]);
    }
  } else {
    MPI_Ssend(v, n*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

int get_vtk_numbodies_mpi (FILE *f)
{
  int N;
  char buf1[512];
  char buf2[512];
  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  if(rank == 0) {
    /* header */
    fgets(buf1,512,f);
    fgets(buf1,512,f);
    fgets(buf1,512,f);

    /* points */
    fgets(buf1,512,f);
    fscanf(f,"%s %d %s",buf1,&N,buf2);
  }

  N /= procs;

  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return N;
}

int read_vtk_file_double_mpi (FILE *f, int n, double x[][3], double v[][3],
                              double m[], double *t, double *dt)
{
  int i,N,n1,n2;
  char buf1[512];
  char buf2[512];
  int rank, procs;
  int p;
  double vec[n][3], sca[n];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  /* header */
  if(rank == 0) {
    fgets(buf1,512,f);
    fscanf(f,"%s %lg %lg\n",buf1,t,dt);
    fgets(buf1,512,f);
  }
  MPI_Bcast(t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* points */
  if(rank == 0) {
    fgets(buf1,512,f);
    fscanf(f,"%s %d %s",buf1,&N,buf2);
  }
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (N != n*procs) {
    printf("data file too large, aborting read\n");
    return -1;
  }
  if(rank == 0) {
    for (i=0; i<n; i++)
      fscanf(f,"%lg %lg %lg\n",&(x[i][0]),&(x[i][1]),&(x[i][2]));
    for (p=1; p < procs; ++p) {
      for (i=0; i<n; i++)
        fscanf(f,"%lg %lg %lg\n",&(vec[i][0]),&(vec[i][1]),&(vec[i][2]));
      MPI_Ssend(vec, n*3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(x, n*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  /* vertices */
  if(rank == 0) {
    fscanf(f,"%s %d %d\n",buf1,&n1,&n2);
    for (i=0; i<n*procs; i++)
      fscanf(f,"%d %d\n",&n1,&n2);
  }

  /* scalar data fields*/
  if(rank == 0) {
    fgets(buf1,512,f);
    fgets(buf1,512,f);
    fgets(buf1,512,f);
    for (i=0; i<n; i++)
      fscanf(f,"%lg\n",&(m[i]));
    for (p=1; p < procs; ++p) {
      for (i=0; i<n; i++)
        fscanf(f,"%lg\n",&(sca[i]));
      MPI_Ssend(sca, n, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(m, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  /* vector data */
  if(rank == 0) {
    fgets(buf1,512,f);
    for (i=0; i<n; i++)
      fscanf(f,"%lg %lg %lg\n",&(v[i][0]),&(v[i][1]),&(v[i][2]));
    for (p=1; p < procs; ++p) {
      for (i=0; i<n; i++)
        fscanf(f,"%lg %lg %lg\n",&(vec[i][0]),&(vec[i][1]),&(vec[i][2]));
      MPI_Ssend(vec, n*3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(v, n*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}
