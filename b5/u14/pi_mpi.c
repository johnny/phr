#include <stdio.h>
#include <math.h>
#include<mpi.h>

double calculate(int rank, int procs, long n){
  double h,sum,x;
  long i,part_length;

  h = 1.0/(double) n;
  sum = 0.0;
  part_length = n/procs;

  if(part_length < 1)
    if(rank==0)
      part_length = n;
    else
      return 0.0;

  for ( i=part_length*rank+1;i<=part_length*(rank+1); i+=1){
    x = h * ((double) i - 0.5);
    sum += (4.0 / (1.0 + x*x));
  }
  return h * sum;
}

/* Berechnung von pi ueber den arcus-tangens */
int main(int argc, char **argv)
{
  int rank, procs;
  const double PI25DT = 3.141592653589793238462643;
  long n, i;
  double pi,pi_part;

  MPI_Status status;
  MPI_Init(&argc, &argv);

  /* Get own rank and total number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  while(1){
    if(rank == 0){
      printf("Anzahl Intervalle (0 beendet): ");
      scanf("%d",&n);

      // teile intervall auf prozesse auf
      for(i = 1; i<procs ; i+=1)
	MPI_Ssend(&n,1,MPI_LONG,i,0,MPI_COMM_WORLD);

    }
    else {
      //empfange anzahl intervalle von prozess 0
      MPI_Recv(&n,1,MPI_LONG,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

    }

    if (n == 0)
      break;
    
    //berechne das teilintervall
    pi_part = calculate(rank,procs,n);

    // sammle ergebnis ein
    MPI_Reduce(&pi_part, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank==0)
      printf("pi: %.16f, Fehler: %.16f\n",
	     pi, fabs(pi - PI25DT));
  }

  MPI_Finalize();

  return 0;
}
