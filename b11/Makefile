CC          = gcc -fopenmp
AR          = ar
AROPTS      = ru
CCOPTS      = -O3 -ffast-math -funroll-loops -fexpensive-optimizations
#CCOPTS      = -g -O0

vanilla_MAIN = nbody_vanilla
vanilla_SOURCES = io_vanilla.o generate_vanilla.o stopwatch.o
vanilla_LIOPTS = -lm

tiled_MAIN = nbody_tiled
tiled_SOURCES = io_vanilla.o generate_vanilla.o stopwatch.o
tiled_LIOPTS = -lm

openmp_MAIN = nbody_openmp
openmp_SOURCES = io_vanilla.o generate_vanilla.o stopwatch.o
openmp_LIOPTS = -lm -lgomp

mpi_MAIN = nbody_mpi
mpi_SOURCES = io_mpi.c generate_mpi.c stopwatch.o
mpi_LIOPTS = -lm -lgomp

all : $(vanilla_MAIN) $(tiled_MAIN) $(openmp_MAIN) $(mpi_MAIN)

$(vanilla_MAIN): $(vanilla_MAIN).o $(vanilla_SOURCES) Makefile
	$(CC) -o $(vanilla_MAIN) $(vanilla_MAIN).o $(vanilla_SOURCES) $(vanilla_LIOPTS)

$(tiled_MAIN): $(tiled_MAIN).o $(tiled_SOURCES) Makefile
	$(CC) -o $(tiled_MAIN) $(tiled_MAIN).o $(tiled_SOURCES) $(tiled_LIOPTS)

$(openmp_MAIN): $(openmp_MAIN).o $(openmp_SOURCES) Makefile
	$(CC) -o $(openmp_MAIN) $(openmp_MAIN).o $(openmp_SOURCES) $(openmp_LIOPTS)

$(mpi_MAIN): $(mpi_MAIN).c $(mpi_SOURCES) Makefile
	mpicc $(CCOPTS) -o $(mpi_MAIN) $(mpi_MAIN).c $(mpi_SOURCES) $(mpi_LIOPTS)

.c.o:
	$(CC) $(CCOPTS) -c $<

.cc.o:
	$(CC) $(CCOPTS) -c $<

clean:
	rm -f $(vanilla_MAIN) $(tiled_MAIN) $(openmp_MAIN) $(mpi_MAIN) *.o
