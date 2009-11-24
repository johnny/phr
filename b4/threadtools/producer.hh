#include "basicthread.hh"
#include "semaphore.hh"

namespace TT {

  class Producer : public TT::BasicThread
  {
  public:
    Producer(double*,Semaphore*,Semaphore*,int,long int);
    ~Producer();
    void run(void);
  private:
    double* buf;
    Semaphore* full;
    Semaphore* empty;
    int k;
    long int runs;
  };
}
