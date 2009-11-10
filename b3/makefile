all: a c stat13 stat102 dyn13 dyn102
a: a2uebung8.c 
	g++ a2uebung8.c -o a -O3 -funroll-loops -fopenmp

stat13: a2uebung8_stat13.c
	g++ a2uebung8_stat13.c -o stat13 -O3 -funroll-loops -fopenmp

stat102: a2uebung8_stat102.c
	g++ a2uebung8_stat102.c -o stat102 -O3 -funroll-loops -fopenmp

dyn13: a2uebung8_dyn13.c
	g++ a2uebung8_dyn13.c -o dyn13 -O3 -funroll-loops -fopenmp

dyn102: a2uebung8_dyn102.c
	g++ a2uebung8_dyn102.c -o dyn102 -O3 -funroll-loops -fopenmp

c: c.c
	g++ c.c -o c -O3 -funroll-loops -fopenmp -g
