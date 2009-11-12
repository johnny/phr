#ifndef __TT_SEMAPHORE_HH__
#define __TT_SEMAPHORE_HH__

#include<iostream>
#include<list>
#include<pthread.h>
#include<semaphore.h>

#include"condition.hh"

namespace TT {
  /** Implements a semaphore using Pthreads condition variables */
  class Semaphore : private Condition<unsigned int>
  {
  public:
    //! make a semaphore with initial value
    Semaphore (int init);

    //! make a semaphore with initial value 0
    Semaphore ();

    //! proceed if value
    void P ();

    //! release semaphore
    void V ();
  };
}

#endif
