#include<iostream>
#include "timer.h"
#include <vector>
#include <stdlib.h>
#include <math.h>

typedef std::vector<double> dv1D;
typedef std::vector<dv1D> dv2D;
typedef std::vector<dv2D> dv3D;

dv3D create_matrix(int n);

double gauss(int n){
  static double sum;
  static int betas[7][3] = {{0,0,1},{0,1,0},{1,0,0},{0,1,1},{1,0,1},{1,1,0},{1,1,1}};
  dv3D U = create_matrix(n);

  struct timeval timer;
  reset_timer(&timer);

  for(int m = 0; m<n;m++){
    for(int a1 = 0; a1<n;a1++)
      for(int a2 = 0; a2<n;a2++)
        for(int a3 = 0; a3<n;a3++){
          sum = 0;
          for(int i = 0; i<7;i++) // Summe in der klammer
              sum += U[((a1-betas[i][0])+n)%n][((a2-betas[i][1])+n)%n][((a3-betas[i][2])+n)%n] + 
		U[((a1+betas[i][0])+n)%n][((a2+betas[i][1])+n)%n][((a3+betas[i][2])+n)%n];

          U[a1][a2][a3] = sum/26;
        }
  }
  return get_timer(timer);
}

double gauss_blocked(int n, int blocksize){
  static double sum;
  static int betas[7][3] = {{0,0,1},{0,1,0},{1,0,0},{0,1,1},{1,0,1},{1,1,0},{1,1,1}};
  dv3D U = create_matrix(n);
  
  struct timeval timer;
  reset_timer(&timer);
for(int a=0;a<n;a+=blocksize){
 for(int b=0;b<n;b+=blocksize){
  for(int c=0;c<n;c+=blocksize){
    for(int m = 0; m<n;m++){
     for(int a1 = a; a1<(a+blocksize);a1++)
       for(int a2 = b; a2<(b+blocksize);a2++)
        for(int a3 = c; a3<(c+blocksize);a3++){
          sum = 0;
          for(int i = 0; i<7;i++) // Summe in der klammer
              sum += U[((a1-betas[i][0])+blocksize)%blocksize][((a2-betas[i][1])+blocksize)%blocksize][((a3-betas[i][2])+blocksize)%blocksize] + 
		U[((a1+betas[i][0])+blocksize)%blocksize][((a2+betas[i][1])+blocksize)%blocksize][((a3+betas[i][2])+blocksize)%blocksize];

          U[a1][a2][a3] = sum/26;
        }
      }
    }
  }
}
  return get_timer(timer);
}

dv3D create_matrix(int n){
  dv3D A(n,dv2D(n,dv1D(n)));

  for(int i = 0; i<n;i++)
    for(int j = 0; j<n;j++)
      for(int k = 0; k<n;k++)
        A[i][j][k]= i*j*k;

  return A;
}

int main(int argc, char* argv[])
{
  double gauss_total_time, gauss_blocked_total_time, gauss_relative_time, gauss_blocked_relative_time;

  long ops;

  int start = 10;//atoi(argv[1]);
  int end = 200;//atoi(argv[2]);
  int runs = 3;//atoi(argv[3]);
  int blocksize = 10;
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
