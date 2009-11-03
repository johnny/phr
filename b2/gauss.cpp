#include<iostream>
#include "timer.h"
#include <stdlib.h>

typedef double*** d3D;
const int betas[7][3] = {{0,0,1},{0,1,0},{1,0,0},{0,1,1},{1,0,1},{1,1,0},{1,1,1}};

d3D create_matrix(int n);

double gauss(int n){
  double sum;
  d3D U = create_matrix(n+2);

  struct timeval timer;
  reset_timer(&timer);

  for(int m = 0; m<n;m++){
    for(int a1 = 1; a1<n+1;a1++)
      for(int a2 = 1; a2<n+1;a2++)
        for(int a3 = 1; a3<n+1;a3++){
          sum = 0;
          for(int i = 0; i<7;i++) // Summe in der klammer
              sum += U[a1-betas[i][0]][a2-betas[i][1]][a3-betas[i][2]] + 
		U[a1+betas[i][0]][a2+betas[i][1]][a3+betas[i][2]];

          U[a1][a2][a3] = sum/26;
        }
  }
  return get_timer(timer);
}

double gauss_blocked(int n, int blocksize){
  double sum;
  d3D U = create_matrix(n+2);
  int a,b,c,m,a1,a2,a3,i;
  
  struct timeval timer;
  reset_timer(&timer);
  for(a=1;a<n+1;a+=blocksize){
    for(b=1;b<n+1;b+=blocksize){
      for(c=1;c<n+1;c+=blocksize){
	for(m = 0; m<n;m++){
	  for(a1 = a; a1<(a+blocksize);a1++)
	    for(a2 = b; a2<(b+blocksize);a2++)
	      for(a3 = c; a3<(c+blocksize);a3++){
		sum = 0;
		for(i = 0; i<7;i++) // Summe in der klammer
		  sum += U[a1-betas[i][0]][a2-betas[i][1]][a3-betas[i][2]] + 
		    U[a1+betas[i][0]][a2+betas[i][1]][a3+betas[i][2]];

		U[a1][a2][a3] = sum/26;
	      }
	}
      }
    }
  }
  return get_timer(timer);
}

d3D create_matrix(int n){
  d3D A;
  int i, j,k;

  A = (double ***)malloc(sizeof(double **) * n);

  for (i = 0 ;  i < n; i++) {
    A[i] = (double **)malloc(sizeof(double *) * n);

    for (j = 0; j < n; j++)
      A[i][j] = (double *)malloc(sizeof(double) * n);
  }

  for(i = 0; i<n;i++)
    for(j = 0; j<n;j++)
      for(k = 0; k<n;k++)
        A[i][j][k]= i*j*k;

  return A;
}

int main(int argc, char* argv[])
{
  double gauss_total_time, gauss_blocked_total_time, gauss_relative_time, gauss_blocked_relative_time;
  int start,end,runs,blocksize;

  long ops;

  if(argc == 5){
    start = atoi(argv[1]);
    end = atoi(argv[2]);
    runs = atoi(argv[3]);
    blocksize = atoi(argv[4]);
  }
  else{
    start = 8;
    end = 120;
    runs = 3;
    blocksize = start;
  }
  for(long n = start; n<end;n+=blocksize){
    gauss_total_time = 0;
    gauss_blocked_total_time = 0;

    for(int i = 0; i<runs;i++){
      gauss_total_time += gauss(n);
      gauss_blocked_total_time += gauss_blocked(n,blocksize);
    }
    gauss_relative_time = gauss_total_time/runs;
    gauss_blocked_relative_time = gauss_blocked_total_time/runs;
    ops = n*n*n*15*n;

    if(false){
      printf("N: %i \n",n);
      printf("Gauss: Total elapsed time: %f s\t",gauss_total_time);
      printf("per run: %f s\t", gauss_relative_time);
      printf("MFlops per run: %f\n", ops/gauss_relative_time/1000000.0);
      printf("Gauss Blocked: Total elapsed time: %f s\t",gauss_blocked_total_time);
      printf("per run: %f s\t", gauss_blocked_relative_time);
      printf("MFlops per run: %f\n", ops/gauss_blocked_relative_time/1000000.0);
    }
    else
      printf("%i %f %f\n", n, ops/gauss_relative_time/1000000.0, ops/gauss_blocked_relative_time/1000000.0);
  }
}
