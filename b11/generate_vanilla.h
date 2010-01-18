#ifndef GENERATE_VANILLA_H
#define GENERATE_VANILLA_H

void cube (int n, long int seed, double size, double m0, double mdelta,  
	   double x[][3], double v[][3], double m[]);
int generate_cube (int n, long int seed, char* name);
void plummer (int n, long int seed,
	      double x[][3], double v[][3], double m[]);
int generate_plummer (int n, long int seed, char* name);

#endif
