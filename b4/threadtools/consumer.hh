#include "basicthread.hh"
#include "semaphore.hh"

namespace TT {

  class Consumer : public TT::BasicThread
  {
    public:
      Consumer(double*,Semaphore*,Semaphore*,int*,int*, int);
      ~Consumer();
      void run();
    private:
      double* buf;
      Semaphore* full;
      Semaphore* empty;
      int* front;
      int* rear;
      int k;
      void manageOrder(double*);
  };
}
