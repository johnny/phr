#include "basicthread.hh"
#include "semaphore.hh"

namespace TT {

  class Producer : public TT::BasicThread
  {
    public:
      Producer(double*,Semaphore*,Semaphore*,int*,int*,int);
      ~Producer();
      void run(void);
      void createOrder(double*);
    private:
      double* buf;
      Semaphore* full;
      Semaphore* empty;
      int* front;
      int* rear;
      int k;
  };
}
