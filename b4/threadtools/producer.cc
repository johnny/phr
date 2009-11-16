#include "producer.hh"

namespace TT {

  Producer::Producer(double* buf,Semaphore* empty,Semaphore* full,int k,long int runs) : BasicThread()
  {
    this->buf = buf;
    this->full = full;
    this->empty = empty;
    this->k = k;
    this->runs = runs;
  }

  Producer::~Producer()
  {

  }

  void Producer::run()
  {
    double i;
    double* position = buf;
    for(i=0;i<runs;++i,++position){
      empty->P();

      if(position == &buf[k]) // Pointer zuruecksetzen
	position = buf;

      *position = i;
      printf("Order sent: %f\n",i);

      full->V();
    }
  }

}
