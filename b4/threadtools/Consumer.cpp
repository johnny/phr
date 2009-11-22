#include "basicthread.hh"
#include "semaphore.hh"

const int N=2;
bool R[N];
TT::Semaphore* res = new TT::Semaphore();//Semaphore der Resourcen
TT::Semaphore* zgr = new TT::Semaphore();//bin.Semaphore sichert Zugriff auf Resourcen ab

class Consumer : public TT::BasicThread
{

public:
	Consumer(int i) : rank(i){
		printf("Consumer created with rank %d\n",i);
	}

	virtual void run() {
       		//Warte auf Resourcen;
		res->P();
       		//Sichere Zugriff auf Resource ab ;
		zgr->P();
       		//Waehle freie Resource i ;
		R[rank]=true;
       		//Entsichere Zugriff ;
		zgr->V();
       		//Arbeite auf Resource i ;
		printf("Resource wird bearbeitet...\n");
       		//Gib Resource i frei ;
		R[rank]=false;
		res->V();
	}

private:
	int rank;
};

int main(int argc, char* argv[])
{
	Consumer* c[N];
	for(int x=0;x<N;x++){
		R[x] = false;
	}
	for(int x=0;x<N;x++){
        	c[x] = new Consumer(x);
		c[x]->start();
	}
	return 0;
}

