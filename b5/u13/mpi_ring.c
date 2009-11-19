#include<mpi.h>
#include<stdio.h>

/* shifting in a ring with colored edges */
int main(int argc, char** argv)
{
  int send, receive;
  int rank, procs;

  int dest, source;

  int tag=50;

  char message[100];

  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /* fill receive buffer */
  receive = rank;

  if (rank!=0)
    {
      sprintf(message,"I am process %d\n",rank);
      dest = 0;
      MPI_Send(message,strlen(message)+1,MPI_CHAR,
	       dest,tag,MPI_COMM_WORLD);
    }
  else
    {
      puts("I am process 0\n");
      for (source=1; source<procs; source++)
	{
	  MPI_Recv(message,100,MPI_CHAR,source,tag,
		   MPI_COMM_WORLD,&status);
	  puts(message);
	}
    }

  /* while(1) */
  /* { */
  /*   send = receive; */

  /*   /\* Do communication *\/ */
  /*   /\* ... *\/ */

  /*   if(rank == 0) */
  /*     printf("Received token of process %d\n", receive); */

  /*   /\* Stop if we got back our own token *\/ */
  /*   if(receive == rank) */
  /*     break; */
  /* } */

  MPI_Finalize();

  return 0;
}
