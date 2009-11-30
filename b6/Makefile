JACOBI = jacobi
JACOBI_SEQ = jacobi_seq

all: jacobi3.c
	mpicc jacobi3.c -o $(JACOBI)

seq: jacobi.c
	gcc jacobi.c -o $(JACOBI_SEQ)

clean:
	rm -f *.o
	rm -f $(JACOBI) $(JACOBI_SEQ)
