#ifndef __TT_MUTEX_HH__
#define __TT_MUTEX_HH__

#include<pthread.h>

namespace TT {

  //! Implements a standard lock
  class Mutex 
  {
  public:
	//! make a mutex in unlocked state
	Mutex ();
	
	//! make a mutex in unlocked state
	~Mutex ();
	
	//! enter critical section
	void lock ();
	
	//! leave critical section
	void unlock ();
	
  private:
	pthread_mutex_t mutex;
  };

}

#endif
