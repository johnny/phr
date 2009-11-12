#ifndef __MY_BASICTHREAD_HH__
#define __MY_BASICTHREAD_HH__

#include<pthread.h>

#include<cstdlib>
#include<iostream>
#include<list>

#include"semaphore.hh"
#include"mutex.hh"
#include"flag.hh"

namespace TT {

  /** ActiveObject provides objects with their own thread of control
      through subclassing.
  */
  class BasicThread
  {
  public:
    /*! Initializes all objects needed in the ActiveObject class. The
      constructor of the derived class will call the start() method
      as its last command in order to start the execution of the
      active object.
    */
    BasicThread ();            

    /*! Destructor of derived class should call stop() on entry in
      order to wait for the active object to complete. Then
      all the private data of the derived class and the base class
      can be destroyed. 
    */
    ~BasicThread ();

    /*! This is the method defining the code of the active object
      and it is implemented in the derived class. The method starts to
      execute as seperate thread by calling the start() method.
    */
    virtual void run () = 0; 


    /*! This method starts the execution of the active object as 
      a seperate thread by starting the action() method with
      pthread_create(). It is usually called in the constructor
      of the derived class as the last command.
    */
    void start ();

    /*! This method waits with pthread_join for the thread of the
      active object to complete. It is usually called as the first
      command in the destructor of the derived class. 
    */
    void stop ();

	/*! returns true when start() has been called and the thread is running
	 */
	bool running ();

	//! copy constructor; works only when thread is not running yet
	BasicThread (const BasicThread&);

	//! assignment operator; works only when thread is not running yet
	BasicThread& operator= (const BasicThread&);


  /*! Send sends a message to the given object. It blocks
      until the receiver has read the message. It is required that
	  o!=this. This implementation is type safe, i.e. it is checked
      whether the corresponding recv uses the same type T. The copy
      uses the operator= of type T and thus can also copy complicated 
      types (unlike an MPI_Send). 

	  All message passing operations can
      only be invoked by active objects (i.e. objects of classes derived
      from ActiveObject).

	  Note also that message passing operations should only be called
      by an object on itself, i.e. not on another object. The idea is
	  that an active object owns a mailbox where other active objects can
      send a message to and the owning object pulls messages from its mailbox.
   */
  template<class T>
  void send (BasicThread& o, T& t);

  /*! Receives a message from the given object. The method blocks until
      the message has arrived. It is required that o!=this. If blocking
      is note desired guard this method by the rprobe method.
   */
  template<class T>
  void recv (BasicThread& o, T& t);

  /*! Blocks until a message from any object is received. This 
      message must be of type T otherwise a type mismatch error
      is reported. A pointer to the sending object is returned as return value.

	  Note that recv should only be called by an object on itself. It
      is assumed that only one thread waits on message arriving at the
      mailbox of an active object.
   */
  template<class T>
  BasicThread& recvany (T& t);

  /*! Returns true if a receive from the given object would not block, 
      i.e. a message has arrived from the object. 
   */
  bool rprobe (BasicThread& o);


  private:
    //! static method that is run as independent thread
    static void* Run (BasicThread* o);

    // thread information
	bool r;              // true when method start() has been called
    Semaphore sema;      // used to join
    pthread_t self;      // my thread id
    pthread_attr_t attr; // attribute

	// mailbox part
	class Envelope {
	public:
	  BasicThread* sender;
	  Flag flag;
	  Envelope(BasicThread* _sender)
	  {
		sender = _sender;
	  }
	  virtual void dummy (void) // enables RTTI
	  {}
	};
	template<class T>
	class TypedEnvelope : public Envelope {
	public:
	  T* data;
	  TypedEnvelope(BasicThread* _sender, T* _data) : Envelope(_sender)
	  {
		data = _data;
	  }
	  virtual void dummy (void)
	  {}
	};

	Mutex mailbox_lock;           // lock for mailbox
	Flag mailbox_notify;          // notify receiving thread that message has arrived
	std::list<Envelope*> mailbox; // list of pending messages for this object
  };

  // message passing implementation

  // send a message to this object
  template<class T>
  void BasicThread::send (BasicThread& o, T& t)
  {
	// check 
	if (&o==this)
	  {
		std::cout << "send(): sender==receiver not allowed" << std::endl;
		return;
	  }

	// make envelope
	TypedEnvelope<T> env(this,&t);

	// put envelope in message queue of receiving object
	o.mailbox_lock.lock();     // get access to mailbox
	o.mailbox.push_back(&env); // put envelope in list
	o.mailbox_lock.unlock();   // release access to mailbox

	// notify receiving thread that new message has arrived
	o.mailbox_notify.signal();

	// block on the envelopes flag until message has been read
	env.flag.wait();
  }

  // receive a message from given object
  template<class T>
  void BasicThread::recv (BasicThread& o, T& t)
  {
	// check 
	if (&o==this)
	  {
		std::cout << "recv(): sender==receiver not allowed" << std::endl;
		return;
	  }

	// make sure that message queue is read on entry
	mailbox_notify.signal(); // ????

	// wait until message has arrived
	std::list<Envelope*>::iterator i;
	while (true)
	  {
		// wait for message to arrive
		mailbox_notify.wait();

		// scan mailbox
		mailbox_lock.lock();
		for (i=mailbox.begin(); i!=mailbox.end(); ++i)
		  if ((*i)->sender==(&o))
			break; // message found
		if (i!=mailbox.end())
		  break;   // leave while loop, stay in lock

		// not found, wait for next message
		mailbox_lock.unlock();
	  }

	// copy message to receivers variable
	TypedEnvelope<T>* p = dynamic_cast<TypedEnvelope<T>*>(*i); // type safe downcast
	if (p!=0)
	  t = *(p->data); // now thats the copy using assignement operator
	else
	  std::cout << "type mismatch in recv" << std::endl;

	// signal to sender that message has been copied
	(*i)->flag.signal();

	// remove envelope from message list
	mailbox.erase(i);

	// release access to message queue
	mailbox_lock.unlock();  
  }

  // wait for message from any active object
  template<class T>
  BasicThread& BasicThread::recvany (T& t)
  {
	// make sure that message queue is read on entry
	mailbox_notify.signal();

	// wait until message has arrived
	std::list<Envelope*>::iterator i;
	while (true)
	  {
		// wait for message to arrive
		mailbox_notify.wait();

		// take first element of nonempty mailbox
		mailbox_lock.lock();
		i=mailbox.begin();
		if (i!=mailbox.end())
		  break;   // leave while loop, stay in lock

		// not found, wait for next message
		mailbox_lock.unlock();
	  }

	// copy message to receivers variable
	TypedEnvelope<T>* p = dynamic_cast<TypedEnvelope<T>*>(*i); // type safe downcast
	if (p!=0)
	  t = *(p->data); // now thats the copy
	else
	  std::cout << "type mismatch in recv" << std::endl;

	// remember sender
	BasicThread* sender=(*i)->sender;

	// signal to sender that message has been copied
	(*i)->flag.signal();

	// remove envelope from message list
	mailbox.erase(i);

	// release access to message queue
	mailbox_lock.unlock();  

	// return sender
	return *sender;
  }


}

#endif
