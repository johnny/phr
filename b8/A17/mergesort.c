#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
void merge(int l, int r, int* a)
{
	int c, i,j,k;
	int m = (r+l)/2;
	int* b = (int *)malloc(sizeof(int)*(r+1));

	for(i=m+1; i>l; i--) b[i-1]=a[i-1];	//linke  Teildatei in Hilfsarray
        for(j=m; j<r; j++) b[r+m-j]=a[j+1];	//rechte Teildatei umgedreht in Hilfsarray
	for(k=l; k<=r; k++)
        	a[k]=(b[i]<b[j])?b[i++]:b[j--];	//füge sortiert in a ein

}

int* Mergesort(int* myList, int n, int rank, int psize, int para)
{
	MPI_Status status;
	if(n==1){
		return myList;
	} else {
		if(rank<psize-2 && para){
			MPI_Ssend(myList,n/2,MPI_INT,rank+2,0,MPI_COMM_WORLD);
		} else Mergesort(myList,n/2,rank,psize,0);
		Mergesort(&myList[n/2],n/2,rank,psize,0);
		if(rank<psize-2 && para){
			MPI_Recv(myList,n/2,MPI_INT,rank+2,0,MPI_COMM_WORLD,&status);
		}
		merge(0,n-1,myList);
		return myList;
	}
}

int main(int argc, char* argv[])
{
	int rank;
	int c;
	int psize;
	int* myList;
	double starttime, endtime;
	int N = atoi(argv[1]);
	myList = (int*)malloc(sizeof(int)*N);

	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&psize);

	if(rank==0){
		// Liste mit Zufallswerten belegen
		srand(time(0));
		for(c=0;c<N;c++){
			myList[c] = rand()%(N+1);
		}
		starttime = MPI_Wtime();
		MPI_Ssend(myList,N/2,MPI_INT,1,0,MPI_COMM_WORLD);
		Mergesort(&myList[N/2],N/2,rank,psize,1);
		MPI_Recv(myList,N/2,MPI_INT,1,0,MPI_COMM_WORLD,&status);
		merge(0,N-1,myList);
		endtime = MPI_Wtime();
		printf("Ergebnis:\n");
		for(c=0;c<N;c++)
			printf("%d%s",myList[c],c!=N-1 ? "," : "\n");
		printf("Laufzeit: %f\n",endtime-starttime);
	} else {
		int size = N/((rank+2)-(rank%2));
		MPI_Recv(myList,size,MPI_INT,rank==1 ? 0 : rank-2,0,MPI_COMM_WORLD,&status);
		printf("(%d) Datenpaket empfangen. Laufe los...\n",rank);
		Mergesort(myList,size,rank,psize,1);
		MPI_Ssend(myList,size,MPI_INT,rank==1 ? 0 : rank-2,0,MPI_COMM_WORLD);
		printf("(%d) Fertig sortiert. Schicke Ergebnis zurück...\n",rank);
	}


	MPI_Finalize();
}
