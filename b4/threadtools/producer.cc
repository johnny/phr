#include "producer.hh"

namespace TT {

  Producer::Producer(double* buf,Semaphore* full,Semaphore* empty,int* front, int* rear, int k) : BasicThread()
  {
    this->buf = buf;
    this->full = full;
    this->empty = empty;
    this->front = front;
    this->rear = rear;
    this->k = k;
  }

  Producer::~Producer()
  {

  }

  void Producer::run()
  {
    double t;
    while(1){
      createOrder(&t);
      empty->P();
      buf[*front] = t;
      *front = (*front+1)%k;
      full->V();
    }
  }

  void Producer::createOrder(double* order)
  {
    *order = 12.0; // welchen Wert nimmt man da? Zufallswert?
    printf("Order sent: %f\n",*order);
  }

}
