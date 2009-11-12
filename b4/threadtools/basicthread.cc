#include<iostream>
#include"basicthread.hh"

namespace TT {

  // constructor
  BasicThread::BasicThread () : sema(0), r(false)
  {
  }

  // destructor
  BasicThread::~BasicThread ()
  {
	if (r)
	  std::cout << "destroying a running thread!" << std::endl;
  }

  // start action method as thread in joinable state
  void BasicThread::start ()
  {
    // create joinable thread
	r = true;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc = pthread_create(&self,&attr,(void*(*)(void*))&BasicThread::Run,(void *)this);
    if (rc)
      {
		std::cout << "ERROR; return code from pthread_create() is " << rc << std::endl;
		exit(-1);
      }
  }

  // wait for thread to complete
  void BasicThread::stop ()
  {
    sema.P(); // sync with end of action()
    pthread_join(self,NULL);
    pthread_attr_destroy(&attr);
	r = false;
  }

  // this is the method started by the thread
  void* BasicThread::Run (BasicThread* o)
  {
    o->run();
    (*o).sema.V(); // sync with stop()
    pthread_exit(0);
  }

  bool BasicThread::running ()
  {
	return r;
  }

  //! copy constructor; works only when thread is not running yet
  BasicThread::BasicThread (const BasicThread& other) : sema(0)
  {
	if (r)
	  {
		std::cout << "attempt to copy a running thread!" << std::endl;
		return;
	  }
	else
	  {
		self = other.self;
		attr = other.attr;
	  }
  }

  //! assignment operator; works only when thread is not running yet
  BasicThread& BasicThread::operator= (const BasicThread& other)
  {
	if (r)
	  {
		std::cout << "attempt to assign a running thread!" << std::endl;
		return *this;
	  }
	else
	  {
		sema = other.sema;
		self = other.self;
		attr = other.attr;
	  }
	return *this;
  }

  // test if message can be received from given object
  bool BasicThread::rprobe (BasicThread& o)
  {
	mailbox_lock.lock();
	for (std::list<Envelope*>::iterator i=mailbox.begin(); i!=mailbox.end(); ++i)
	  if ((*i)->sender==(&o))
		{
		  mailbox_lock.unlock();
		  return true;
		}
	mailbox_lock.unlock();
	return false;
  }

}
