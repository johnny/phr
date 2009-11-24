#include"semaphore.hh"

namespace TT {

  // implementation of the alternative semaphore
  Semaphore::Semaphore (int init) : Condition<unsigned int>(init)
  {
  }

  Semaphore::Semaphore () : Condition<unsigned int>(0)
  {
  }

  void Semaphore::P ()
  {
    aquire();

    while(!value())
      wait();

    --value();

    release();  

  }

  void Semaphore::V ()
  {
    aquire();

    ++value();

    signal();

    release();
  }

}
