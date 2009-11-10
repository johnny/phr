#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

int main(int argc, char** argv)
{
	int n;
	for(n=1;n<4;n++){
		int N=128*(int)pow(2,n);
		int p;
                int** mat1=(int **) malloc(sizeof(int*)*N);
                int** mat2=(int **) malloc(sizeof(int*)*N);
                int** mat3=(int **) malloc(sizeof(int*)*N);
		for(p=0;p<N;p++){
			mat1[p]=(int *) malloc(sizeof(int)*N);
			mat2[p]=(int *) malloc(sizeof(int)*N);
			mat3[p]=(int *) malloc(sizeof(int)*N);
		}
		int k,i,j,c,m;
		double tall;
		for(m=0;m<11;m++){
			omp_set_num_threads((int)pow(2,m));

			// Matrizen initialisieren
			for(j=0;j<N;j++){
				for(k=0;k<N;k++){
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
				#pragma omp parallel for schedule (static)
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
			printf("m: %d n: %d time: %f\n",(int)pow(2,m),128*pow(2,n),tall/c);
		}
		free(mat1);
		free(mat2);
		free(mat3);
	}
	return 0;
}
