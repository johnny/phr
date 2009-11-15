#include "producer.hh"
#include "consumer.hh"
#include "semaphore.hh"

  int main(int argc, char* argv[])
  {
    int k = 20;
    TT::Semaphore empty = TT::Semaphore(k);
    TT::Semaphore full = TT::Semaphore(0);
    double* buf = (double*)malloc(sizeof(double)*k);
    int front = 0;
    int rear = 0;
    TT::Producer p = TT::Producer(buf,&empty,&full,&front,&rear,k);
    TT::Consumer c = TT::Consumer(buf,&empty,&full,&front,&rear,k);
    p.start();
    c.start();
  }

