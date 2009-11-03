#include<iostream>
#include "timer.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

typedef double*** d3D;
const int betas[7][3] = {{0,0,1},{0,1,0},{1,0,0},{0,1,1},{1,0,1},{1,1,0},{1,1,1}};

d3D create_matrix(int n);

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
  double gauss_blocked_total_time, gauss_blocked_relative_time;
  int start,end,runs,stepsize;

  long ops;
  long n = 120;

  if(argc == 2)
    runs = atoi(argv[1]);
  else
    runs = 3;

  for(int i = 5; i<=5; i*=5)
    for(int j = 3; j<=3; j*=3)
      for(int k = 4; k<=4; k*=2){
	int blocksize = i*j*k;
	gauss_blocked_total_time = 0;

	for(int s = 0; s<runs;s++){
	  gauss_blocked_total_time += gauss_blocked(n,blocksize);
	}
	gauss_blocked_relative_time = gauss_blocked_total_time/runs;
	ops = n*n*n*15*n;
	printf("%i %f\n", blocksize, ops/gauss_blocked_relative_time/1000000.0);
      }
}
