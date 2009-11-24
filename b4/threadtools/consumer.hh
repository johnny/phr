#include "basicthread.hh"
#include "semaphore.hh"

namespace TT {

  class Consumer : public TT::BasicThread
  {
  public:
    Consumer(double*,Semaphore*,Semaphore*,int,long int);
    ~Consumer();
    void run();
  private:
    double* buf;
    Semaphore* full;
    Semaphore* empty;
    int k;
    long int runs;
  };
}
