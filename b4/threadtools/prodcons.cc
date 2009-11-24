#include "producer.hh"
#include "consumer.hh"
#include "semaphore.hh"

  int main(int argc, char* argv[])
  {
    int k = 500;
    long int runs = 100000;
    TT::Semaphore empty = TT::Semaphore(k);
    TT::Semaphore full = TT::Semaphore(0);
    double buf[k];
    TT::Producer p = TT::Producer(buf,&empty,&full,k,runs);
    TT::Consumer c = TT::Consumer(buf,&empty,&full,k,runs);
    c.start();
    p.start();
    c.stop();
    p.stop();
  }

