#ifndef IO_MPI_H
#define IO_MPI_H

void write_vtk_file_double_mpi (FILE *f, int n, double x[][3], double v[][3],
                                double m[], double t, double dt);
int get_vtk_numbodies_mpi (FILE *f);
int read_vtk_file_double_mpi (FILE *f, int n, double x[][3], double v[][3],
                              double m[], double *t, double *dt);

#endif
