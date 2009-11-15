#include "consumer.hh"

namespace TT {

  Consumer::Consumer(double* buf,Semaphore* full,Semaphore* empty,int* front, int* rear, int k) : BasicThread()
  {
    this->buf = buf;
    this->full = full;
    this->empty = empty;
    this->front = front;
    this->rear = rear;
    this->k = k;
  }

  Consumer::~Consumer()
  {

  }

  void Consumer::run()
  {
    double t;
    while(1){
      full->P();
      t = buf[*rear];
      *rear = (*rear+1)%k;
      empty->V();
      manageOrder(&t);
    }
  }

  void Consumer::manageOrder(double* order){
    printf("Order received: %f\n",*order);
  }

}
