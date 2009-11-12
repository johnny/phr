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
    printf("P inside\n");
    aquire();
    printf("P aquired\n");
    printf("\t P: value before %i\n",value());
    while(!value()){
      printf("\t has to wait %i\n",value());
      wait();
      printf("\t should go %i\n",value());
      // aquire();
      printf("\t should go after aquired %i\n",value());
    }
    printf("\t P decrements\n");
    --value();
    printf("\t P: value after %i\n",value());
    release();  
  }

  void Semaphore::V ()
  {
    printf("V inside\n");
    aquire();
    printf("V aquired\n");
    if(!value()){
      printf("\t V signals\n");
      signal();
    }
    printf("\t V increments\n");
    printf("\t value before %i\n",value());
    ++value();
    printf("\t value after %i\n",value());
    release();
  }

}
