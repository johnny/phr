all: jm js
jm: jacobi_mpi.c
	mpicxx -g -o jm jacobi_mpi.c
js: jacobi.c
	mpicxx -o js jacobi.c
