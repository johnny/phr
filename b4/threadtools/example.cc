#include<iostream>
#include"tt.hh"

TT::Mutex lock;
TT::Flag flag; 

// a minimalistic thread
class Minimalistic : public TT::BasicThread
{
public:
  virtual void run () {
	flag.wait();
	lock.lock();
    std::cout << "thread running" << std::endl;
	lock.unlock();
  }
};

// a thread with parameters starting itself automatically
class SelfStartingWithArgs : public TT::BasicThread
{
public:
  SelfStartingWithArgs (int i) : rank(i) {
	this->start();
  }
  ~SelfStartingWithArgs () {
	this->stop();
  }
  virtual void run ()
  {
	lock.lock();
    std::cout << "thread " << rank << " running" << std::endl;
	lock.unlock();
	flag.signal();
  }
private:
  int rank;
};

// element of a pool of threads with argument
class PoolElement : public TT::ThreadPoolElement<PoolElement>
{
public:
  PoolElement (TT::Barrier& b) : barrier(b) {}
  virtual void run ()
  {
	barrier.sync();
    std::cout << "array thread " << context().rank() 
			  << " of " << context().size() << std::endl;
    PoolElement& next = 
	  context().pool()[(context().rank()+1)%context().size()];
  }
private:
  int arg;
  TT::Barrier& barrier;
};


int main (int argc, char *argv[])
{
  Minimalistic a; a.start();
  SelfStartingWithArgs b(5);
  TT::Barrier barrier(10);
  TT::ThreadPool<PoolElement> c(10,PoolElement(barrier));
  c.start();
  a.stop();
  c.stop();
}
