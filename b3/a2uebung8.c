#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#ifdef _OPENMP
  #include <omp.h>
#endif

int main(int argc, char** argv)
{
  int n;
  for(n=0;n<3;n++){
    int N=128*(int)pow(2,n);
    int p;
    double** mat1=(double **) malloc(sizeof(double*)*N);
    double** mat2=(double **) malloc(sizeof(double*)*N);
    double** mat3=(double **) malloc(sizeof(double*)*N);
    for(p=0;p<N;p++){
      mat1[p]=(double *) malloc(sizeof(double)*N);
      mat2[p]=(double *) malloc(sizeof(double)*N);
      mat3[p]=(double *) malloc(sizeof(double)*N);
    }
    int k,i,j,c,m;
    double tall;
    for(m=0;m<11;m++){
      omp_set_num_threads((int)pow(2,m));
      tall = 0.0;

      // Mehrfacher Schleifendurchlauf fuer gemitteltes Ergebnis
      for(c=0;c<5;c++){

	// Matrizen initialisieren
	for(j=0;j<N;j++){
	  for(k=0;k<N;k++){
	    mat1[j][k]=(double)k+j;
	    mat2[j][k]=(double)k+j;
	    mat3[j][k]= 0.0;
	  }
	}
	
	// Timer setzen
	double tstart = omp_get_wtime();
	// Matrizen multiplizieren
        #pragma omp parallel for schedule (static) private(j,k) shared(mat3,mat1,mat2)
	for(i=0;i<N;i++){
	  for(j=0;j<N;j++){
	    for(k=0;k<N;k++){
	      mat3[i][j]+=mat1[i][k]*mat2[k][j];
	    }
	  }
	}

	//Timer stoppen und ausgeben
	double tend = omp_get_wtime();
	tall += tend-tstart;
      }
      if(false)
	printf("m: %d n: %d time: %f\n",(int)pow(2,m),128*(int)pow(2,n),tall/c);
      else
	printf("%d %4d %f\n",(int)pow(2,m),128*(int)pow(2,n),tall/c);
    }
    free(mat1);
    free(mat2);
    free(mat3);
  }
  return 0;
}
