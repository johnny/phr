#ifndef __TT_BARRIER_HH__
#define __TT_BARRIER_HH__

#include"mutex.hh"
#include"semaphore.hh"

namespace TT {

  //! Implements a barrier for a fixed number of processes
  class Barrier 
  {
  public:
	//! make a barrier for a given number of threads
	Barrier (unsigned int nthreads);

	//! perform barrier synchronization
	void sync ();

  private:
	unsigned int n;
	unsigned int turn;
	Mutex mutex;
	unsigned int count;
	Semaphore sem[2];
  };
}

#endif
