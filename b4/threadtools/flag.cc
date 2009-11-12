#include"flag.hh"

namespace TT {

  // Flag
  Flag::Flag () : Condition<bool>(false) 
  {
  }

  Flag::~Flag ()
  {
  }

  void Flag::signal ()
  {
	this->aquire();
	this->value() = true;
	Condition<bool>::signal();
	this->release();
  }

  void Flag::wait ()
  {
	this->aquire();
	while (this->value()==false)
	  Condition<bool>::wait();
	this->value() = false;
	this->release();
  }

}
