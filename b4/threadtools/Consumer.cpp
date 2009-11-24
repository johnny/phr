#include "basicthread.hh"
#include "semaphore.hh"

const int N=22;
bool R[2] = {false,false};
TT::Semaphore* res = new TT::Semaphore(2);//Semaphore der Resourcen
TT::Mutex* zgr = new TT::Mutex();//bin.Semaphore sichert Zugriff auf Resourcen ab

class Consumer : public TT::BasicThread
{

public:
	Consumer(int k) : rank(k){
		printf("Consumer created with rank %d\n",i);
	}

	void run() {
		while(1){
	       		//Warte auf Resourcen;
			res->P();
       			//Sichere Zugriff auf Resource ab ;
			zgr->lock();
       			//Waehle freie Resource i ;
			if(R[0]){
				if(R[1]){
					printf("Beide Resourcen belegt, Pech gehabt!\n");
				} else {
					R[1]=true;
					i=1;
				}
			} else {
				R[0]=true;
				i=0;
			}
	       		//Entsichere Zugriff ;
			zgr->unlock();
       			//Arbeite auf Resource i ;
			printf("Thread[%d]: Resource[%d] wird bearbeitet...\n",rank,i);
       			//Gib Resource i frei ;
			R[i]=false;
			res->V();
		}
	}

private:
	int rank;
	int i;
};

int main(int argc, char* argv[])
{
	Consumer* c[N];
	for(int x=0;x<N;x++){
        	c[x] = new Consumer(x);
		c[x]->start();
	}
	while(1); // Absichern, dass Elternthread nicht Kinderthreads beendet
	return 0;
}

