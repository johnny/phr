#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

int main(int argc, char** argv)
{
		int N=500;
		int p;
                int** mat1=(int **) malloc(sizeof(int*)*N);
                int** mat2=(int **) malloc(sizeof(int*)*N);
                int** mat3=(int **) malloc(sizeof(int*)*N);
		for(p=0;p<N;p++){
			mat1[p]=(int *) malloc(sizeof(int)*N);
			mat2[p]=(int *) malloc(sizeof(int)*N);
			mat3[p]=(int *) malloc(sizeof(int)*N);
		}
		int c,m;
		double tall;
		for(m=1;m<=16;m*=2){
			omp_set_num_threads(m);

			// Matrizen initialisieren
			for(int j=0;j<N;j++){
				for(int k=0;k<N;k++){
					mat1[j][k]=k+j;
					mat2[j][k]=k+j;
					mat3[j][k]=k+j;
				}
			}

			// Mehrfacher Schleifendurchlauf fuer gemitteltes Ergebnis
			for(c=0;c<5;c++){

				// Timer setzen
				double tstart = omp_get_wtime();
				// Matrizen multiplizieren
				#pragma omp parallel for schedule (static,102)
				for(int i=0;i<N;i++){
					for(int j=0;j<N;j++){
						for(int k=0;k<N;k++){
							mat3[i][j]+=mat1[i][k]*mat2[k][j];
						}
					}
				}

				//Timer stoppen und ausgeben
				double tend = omp_get_wtime();
				tall += tend-tstart;
			}
			printf("static,102 - m: %d n: %d time: %f\n",m,500,tall/c);
		}
		free(mat1);
		free(mat2);
		free(mat3);
	return 0;
}
