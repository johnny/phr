#ifndef __TT_FLAG_HH__
#define __TT_FLAG_HH__

#include"condition.hh"

namespace TT {

  //! Implements a flag that can be used to wait for arrival of another thread
  class Flag : Condition<bool>
  {
  public:
	//! set up flag variable
	Flag();

	//! signal to other process that on has arrived
	void signal ();

	//! wait for someone to arrive
	void wait ();

	//! destroy flag variable
	~Flag();
  };

}

#endif
