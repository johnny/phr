#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif

int main(int argc, char** argv)
{
  int n,p;
  const int max=10;
  double tall[max];
  const int number_of_threads = 4;
  omp_set_num_threads(number_of_threads);
  double** mat1,** mat2,** mat3;
  double tstart, tend;
  int k,i,j,c,N;
  
#pragma omp parallel for schedule (static) private (N,mat1,mat2,mat3,j,k,c,i,n,tstart,tend)
  for(n = 1; n<max;n++){
    N=(int)pow(2,n);

    mat1=(double **) malloc(sizeof(double*)*N);
    mat2=(double **) malloc(sizeof(double*)*N);
    mat3=(double **) malloc(sizeof(double*)*N);

    for(p=0;p<N;p++){
      mat1[p]=(double *) malloc(sizeof(double)*N);
      mat2[p]=(double *) malloc(sizeof(double)*N);
      mat3[p]=(double *) malloc(sizeof(double)*N);
    }

    for(c=0;c<5;c++){
      tall[n] = 0.0;

      // Matrizen initialisieren
      for(j=0;j<N;j++){
	for(k=0;k<N;k++){
	  mat1[j][k]=(double)k+j;
	  mat2[j][k]=(double)k+j;
	  mat3[j][k]= 0.0;
	}
      }
	
      // Timer setzen
      tstart = omp_get_wtime();
      // Matrizen multiplizieren
      for(i=0;i<N;i++){
	for(j=0;j<N;j++){
	  for(k=0;k<N;k++){
	    mat3[i][j]+=mat1[i][k]*mat2[k][j];
	  }
	}
      }

      //Timer stoppen und ausgeben
      tend = omp_get_wtime();
      tall[n] += tend-tstart;
    }
    if(false)
      printf("n: %d time: %f\n",(int)pow(2,n),tall[n]/c);
    else
      printf("%4d %f\n",(int)pow(2,n),tall[n]/c);
    free(mat1);
    free(mat2);
    free(mat3);
  }
  return 0;
}
