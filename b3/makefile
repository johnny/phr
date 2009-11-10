all: a b test
a: a2uebung8.c 
	g++ a2uebung8.c -o a -O3 -funroll-loops -fopenmp

b: a2uebung8_dyn.c
	g++ a2uebung8_dyn.c -o b -O3 -funroll-loops -fopenmp

test: example.c
	g++ example.c -o test -O3 -funroll-loops -fopenmp