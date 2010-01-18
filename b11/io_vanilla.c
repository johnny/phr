#include <stdio.h>


void write_vtk_file_double (FILE *f, int n, double x[][3], double v[][3], double m[],
			    double t, double dt)
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
    fprintf(f,"%22.16g %22.16g %22.16g\n",x[i][0],x[i][1],x[i][2]);

  /* vertices */
  fprintf(f,"%s %d %d\n","VERTICES",n,2*n);
  for (i=0; i<n; i++)
    fprintf(f,"%d %d\n",1,i);

  /* scalar data fields*/
  fprintf(f,"%s %d\n","POINT_DATA",n);
  fprintf(f,"%s\n","SCALARS mass float");
  fprintf(f,"%s\n","LOOKUP_TABLE default");
  for (i=0; i<n; i++)
    fprintf(f,"%22.16g\n",m[i]);

  /* vector data */
  fprintf(f,"%s\n","VECTORS velocity float");
  for (i=0; i<n; i++)
    fprintf(f,"%22.16g %22.16g %22.16g\n",v[i][0],v[i][1],v[i][2]);
}

int get_vtk_numbodies (FILE *f)
{
  int N;
  char buf1[512];
  char buf2[512];

  /* header */
  fgets(buf1,512,f);
  fgets(buf1,512,f);
  fgets(buf1,512,f);

  /* points */
  fgets(buf1,512,f);
  fscanf(f,"%s %d %s",buf1,&N,buf2);
  return N;
}

int read_vtk_file_double (FILE *f, int n, double x[][3], double v[][3], double m[],
			  double *t, double *dt)
{
  int i,N,n1,n2;
  char buf1[512];
  char buf2[512];

  /* header */
  fgets(buf1,512,f);
  fscanf(f,"%s %lg %lg\n",buf1,t,dt);
  fgets(buf1,512,f);

  /* points */
  fgets(buf1,512,f);
  fscanf(f,"%s %d %s",buf1,&N,buf2);
  if (N>n) 
    {
      printf("data file too large, aborting read\n");
      return -1;
    }
  n = N;
  for (i=0; i<n; i++)
    fscanf(f,"%lg %lg %lg\n",&(x[i][0]),&(x[i][1]),&(x[i][2]));

  /* vertices */
  fscanf(f,"%s %d %d\n",buf1,&n1,&n2);
  for (i=0; i<n; i++)
    fscanf(f,"%d %d\n",&n1,&n2);

  /* scalar data fields*/
  fgets(buf1,512,f);
  fgets(buf1,512,f);
  fgets(buf1,512,f);
  for (i=0; i<n; i++)
    fscanf(f,"%lg\n",&(m[i]));

  /* vector data */
  fgets(buf1,512,f);
  for (i=0; i<n; i++)
    fscanf(f,"%lg %lg %lg\n",&(v[i][0]),&(v[i][1]),&(v[i][2]));
}
