#ifndef __MY_THREADPOOL_HH__
#define __MY_THREADPOOL_HH__

#include<iostream>
#include<vector>

#include"basicthread.hh"

namespace TT {

  template<class T> class ThreadPoolElement;
  template<class T> class ThreadPoolContext;
  template<class T> class ThreadPool;

  template<class T>
  class ThreadPoolElement : public BasicThread
  {
  public:
	typedef T ThreadPoolElementType;
	typedef ThreadPool<T> ThreadPoolType;
    friend class ThreadPool<T>;

	ThreadPoolContext<T>& context ()
	{
	  return context_;
	}

  private:
    ThreadPoolContext<T> context_;
  };

  template<class T>
  class ThreadPoolContext 
  {
  public:
	typedef T ThreadPoolElementType;
	typedef ThreadPool<T> ThreadPoolType;
    friend class ThreadPool<T>;

    ThreadPoolElementType& peer (int i) const
    {
      return (*p)[i];
    }
    int rank () const
    {
      return i;
    }
    int size () const
    {
      return p->size();
    }
    ThreadPoolType& pool ()
    {
      return *p;
    }
  private:
    int i;
    ThreadPoolType *p;
  };


  template<class T>
  class ThreadPool
  {
  public:
    ThreadPool (int P, const T& t) : pool(P,t)
    {
      for (int i=0; i<P; i++)
		{
		  pool[i].context_.i = i;
		  pool[i].context_.p = this;
		}
    }

    int size () const
    {
      return pool.size();
    }

    void start ()
    {
      for (unsigned int i=0; i<pool.size(); i++)
		pool[i].start();
    }

    void stop ()
    {
      for (unsigned int i=0; i<pool.size(); i++)
		pool[i].stop();
    }

    T& operator[] (int i)
    {
      return pool[i];
    }

    const T& operator[] (int i) const
    {
      return pool[i];
    }

  private:
	std::vector<T> pool;
    ThreadPool (const ThreadPool&) {}
    ThreadPool& operator= (const ThreadPool&) {}
  };
}

#endif
