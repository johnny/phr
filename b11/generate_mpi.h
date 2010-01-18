#ifndef GENERATE_MPI_H
#define GENERATE_MPI_H

void cube_mpi (int n, long int seed, double size, double m0, double mdelta,  
               double x[][3], double v[][3], double m[]);
void plummer_mpi (int n, long int seed,
                  double x[][3], double v[][3], double m[]);

#endif
