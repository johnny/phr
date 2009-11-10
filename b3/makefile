all: a b c d e
a: a2uebung8.c 
	g++ a2uebung8.c -o a -O3 -funroll-loops -fopenmp

b: a2uebung8_stat13.c
	g++ a2uebung8_stat13.c -o stat13 -O3 -funroll-loops -fopenmp

c: a2uebung8_stat102.c
	g++ a2uebung8_stat102.c -o stat102 -O3 -funroll-loops -fopenmp

d: a2uebung8_dyn13.c
	g++ a2uebung8_dyn13.c -o dyn13 -O3 -funroll-loops -fopenmp

e: a2uebung8_dyn102.c
	g++ a2uebung8_dyn102.c -o dyn102 -O3 -funroll-loops -fopenmp

test: example.c
	g++ example.c -o test -O3 -funroll-loops -fopenmp
