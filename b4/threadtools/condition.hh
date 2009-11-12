#ifndef __MY_CONDITION_HH__
#define __MY_CONDITION_HH__

#include<pthread.h>
#include<semaphore.h>

#include<iostream>
#include<list>

namespace TT {
  /** A class for condition synchronisation. A Condition combines
      a variable of the condent type T, a mutex variable and a condition
      variable
  */
  template<class T>
  class Condition
  {
  public:
    //! constructor with initial value for data
    Condition (const T& x);

    //! constructor with default values (no initialization for built-in types)
    Condition ();

    //! aquire exclusive access to the lock
    void aquire ();

    //! access data to check or modify it
    T& value ();

    //! access data to check or modify it
    const T& value () const;

    //! unblock a waiting thread and KEEP the lock
    void signal ();

    //! unblock all waiting threads and KEEP the lock
    void broadcast ();

    //! wait for signal and release the lock (atomically)
    void wait ();

    //! release exclusive access to the lock
    void release ();

    //! release all resources
    ~Condition ();
	
  private:
    T data;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    unsigned int waiting;
  };

  template<class T>
  Condition<T>::Condition (const T& x)
  {
    data = x;
    pthread_mutex_init(&mutex,NULL);
    pthread_cond_init(&cond,NULL);
    waiting = 0;
  }

  template<class T>
  Condition<T>::Condition ()
  {
    pthread_mutex_init(&mutex,NULL);
    pthread_cond_init(&cond,NULL);  
    waiting = 0;
  }

  template<class T>
  Condition<T>::~Condition ()
  {
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);  
  }

  template<class T>
  void Condition<T>::aquire ()
  {
    pthread_mutex_lock(&mutex);
  }

  template<class T>
  void Condition<T>::release ()
  {
    pthread_mutex_unlock(&mutex);
  }

  template<class T>
  T& Condition<T>::value ()
  {
    return data;
  }

  template<class T>
  const T& Condition<T>::value () const
  {
    return data;
  }

  template<class T>
  void Condition<T>::signal ()
  {
    if (waiting>0)
      {
		waiting--;
		pthread_cond_signal(&cond);
      }
  }

  template<class T>
  void Condition<T>::broadcast ()
  {
    if (waiting>0)
      {
		waiting=0;
		pthread_cond_broadcast(&cond);
      }
  }

  template<class T>
  void Condition<T>::wait ()
  {
    waiting++;
    pthread_cond_wait(&cond,&mutex);
  }
}

#endif
