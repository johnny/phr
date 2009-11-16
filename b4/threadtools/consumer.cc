#include "consumer.hh"

namespace TT {

  Consumer::Consumer(double* buf,Semaphore* empty,Semaphore* full,int k,long int runs) : BasicThread()
  {
    this->buf = buf;
    this->full = full;
    this->empty = empty;
    this->k = k;
    this->runs = runs;
  }

  Consumer::~Consumer()
  {

  }

  void Consumer::run()
  {
    double* position = buf;
    for(double i=0;i<runs;++i,++position){
      full->P();

      if(position == &buf[k]) // Pointer zuruecksetzen
	position = buf;

      printf("Order received: %f\n",*position);

      empty->V();
    }
  }

}
