#include<iostream>
#include "timer.h"
#include <vector>
#include <stdlib.h>

typedef std::vector<double> dv1D;
typedef std::vector<dv1D> dv2D;
typedef std::vector<dv2D> dv3D;

double gauss(int n){
  static double sum;
  static int betas[7][3] = {{0,0,1},{0,1,0},{1,0,0},{0,1,1},{1,0,1},{1,1,0},{1,1,1}};
  int misses;  // nur temporaer drin. ich will versuchen direkt die Anzahl der Operationen anzugeben. Siehe ops in main
  dv3D U = create_matrix(n);

  struct timeval timer;
  reset_timer(&timer);

  for(int m = 0; m<n;m++){
    misses = 0;
    for(int a1 = 0; a1<n;a1++)
      for(int a2 = 0; a2<n;a2++)
        for(int a3 = 0; a3<n;a3++){
          sum = 0;
          for(int i = 0; i<7;i++) // Summe in der klammer
            if( a1-betas[i][0] >= 0 && a2-betas[i][1] >= 0 && a3-betas[i][2] >= 0 && // sonst gibt es segfaults
                a1+betas[i][0] < n && a2+betas[i][1] < n && a3+betas[i][2] < n)
              sum += U[a1-betas[i][0]][a2-betas[i][1]][a3-betas[i][2]] + U[a1+betas[i][0]][a2+betas[i][1]][a3+betas[i][2]];
            else
              misses++; // zaehle die nicht beruecksichtigten terme mit

          U[a1][a2][a3] = sum/26;
        }
  }
  std::cout << misses << std::endl;
  return get_timer(timer);
}

double gauss_blocked(dv3D U){
  return 0; // TODO
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
  int end = 41;//atoi(argv[2]);
  int runs = 1;//atoi(argv[3]);
  int stepwidth = 10;

  for(long n = start; n<end;n+=stepwidth){
    gauss_total_time = 0;
    gauss_blocked_total_time = 0;

    for(int i = 0; i<runs;i++){
      gauss_total_time += gauss(n);
      gauss_blocked_total_time += gauss_blocked(n);
    }
    gauss_relative_time = gauss_total_time/runs;
    gauss_blocked_relative_time = gauss_blocked_total_time/runs;
    ops = 2*n*n*n*n;

    if(true){
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
