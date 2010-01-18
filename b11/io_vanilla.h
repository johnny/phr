#ifndef IO_VANILLA_H
#define IO_VANILLA_H

void write_vtk_file_double (FILE *f, int n, double x[][3], double v[][3], double m[], 
			    double t, double dt);
int get_vtk_numbodies (FILE *f);
int read_vtk_file_double (FILE *f, int n, double x[][3], double v[][3], double m[], 
			  double *t, double *dt);

#endif
