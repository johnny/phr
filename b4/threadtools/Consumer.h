#include "basicthread.hh"
#include "semaphore.hh"

const int N=2;
bool R[N];
TT::Semaphore* res = new TT::Semaphore();//Semaphore der Resourcen
TT::Semaphore* zgr = new TT::Semaphore();//bin.Semaphore sichert Zugriff auf Resourcen ab

class Consumer : public TT::BasicThread
{
        private:
                static int pID;
	public:
		Consumer(int i) : rank(i);
		~Consumer();
		void run();
}
