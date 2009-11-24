all: pim pis gnuplot
pim: pi_mpi.c
	mpicc -o pim pi_mpi.c
pis: piseq.c
	mpicc -o pis piseq.c
gnuplot: plot.dat
	gnuplot -persist 'plot'
