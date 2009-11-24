#include<mpi.h>
#include<stdio.h>
/* shifting in a ring with colored edges */
int main(int argc, char** argv)
{
  int receive[2];
  int send, breaker;
  int rank, procs;
  int startflag=1; // signals process with even rank
  MPI_Status status;
  MPI_Init(&argc, &argv);

  /* Get own rank and total number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* fill receive buffer */
  receive[1] = rank;
  receive[0] = rank;
  while(1)
  {
        //Send
        if(rank%2==0)startflag=0;
        if(startflag==0){
                printf("%d: Sending...%d\n",rank,receive[1]);
                if(rank==procs-1){
                        MPI_Ssend(&receive[1],1,MPI_INT,0,0,MPI_COMM_WORLD);
                } else {
                        MPI_Ssend(&receive[1],1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
                }
        }
        if(breaker==1)break;
        if(rank%2==1)receive[1]=receive[0];
        //Receive
        if(rank!=0){
                MPI_Recv(rank%2==0 ? &receive[1] : &receive[0],1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        } else {
                MPI_Recv(&receive[1],1,MPI_INT,procs-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        }
        startflag=0;

      printf("%d: Received token %d\n", rank, receive[1]);

    /* Stop if we got back our own token */
    if((rank%2==0 && receive[1] == rank) || (rank%2==1 && receive[0] == rank)){
      printf("%d: bin fertig!\n",rank);
      if(rank%2==0)break;
        else breaker=1;
    }
  }

  MPI_Finalize();

  return 0;
}
