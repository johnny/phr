#include "barrier.hh"

namespace TT {

  // Barrier implementation
  Barrier::Barrier (unsigned int nthreads)
  {
	n = nthreads;
	turn = 0;
	count = 0;
  }

  void Barrier::sync ()
  {
	mutex.lock();          // enter lock
	count++;               // increment counter for current turn
	if (count<n)           // check if not the last one
	  {
		int myturn=turn;   // read turn before unlocking
		mutex.unlock();    // leave lock
		sem[myturn].P();   // wait on semaphore for this turn
	  }
	else
	  {
		for (int i=0; i<n-1; i++)
		  sem[turn].V();   // wake up waiting n-1 threads
		count = 0;         // reset counter, we are still in the lock!
		turn = 1-turn;     // change to other semaphore
		mutex.unlock();    // leave lock, now others may enter
	  }
  }

}
