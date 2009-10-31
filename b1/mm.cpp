#include<iostream>
#include "timer.h"
#include <vector>
#include <stdlib.h>

using namespace std;

double run_calculation(int n, vector<vector<double> >& A, vector<vector<double> >& B){
  vector<vector<double> > C(n,vector<double>(n));
  struct timeval timer;
  reset_timer(&timer);
  
  for(int i = 0; i<n;i++){
    for(int j = 0; j<n;j++){
      for(int k = 0;k<n;k++){
	C[i][j] = A[i][k]*B[k][j];
      }
    }
  }

  return get_timer(timer);
}

vector<vector<double> > create_matrix(int n){
  vector<vector<double> > A(n,vector<double>(n));
  
  for(int i = 0; i<n;i++){
    for(int j = 0; j<n;j++){
      A[i][j]= i*j;
    }
  }
  return A;
}

int main(int argc, char* argv[])
{
  double total_time;
  double relative_time;
  long ops;
  
  int start = atoi(argv[1]);
  int end = atoi(argv[2]);
  int runs = atoi(argv[3]);
  int stepwidth = 2;

  for(long n = start; n<end;n+=stepwidth){
    total_time = 0;

    vector<vector<double> > A = create_matrix(n);
    vector<vector<double> > B = create_matrix(n);

    for(int i = 0; i<runs;i++){
      total_time += run_calculation(n,A,B);
    }
    relative_time = total_time/runs;
    ops = 2*n*n*n;
    
    if(false){
      printf("N: %i \t",n);
      printf("Total elapsed time: %f s\t",total_time);
      printf("per run: %f s\t", relative_time);
      printf("MFlops per run: %f\n", ops/relative_time/1000000.0);
    }
    else
      printf("%i %f\n", n, ops/relative_time/1000000.0);
  }
}

