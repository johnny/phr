#include<iostream>
#include"tt.hh"

// element of a pool of threads with argument
class MyThread : public TT::ThreadPoolElement<MyThread>
{
public:
  virtual void run ()
  {
	int rank=context().rank();
	int size=context().size();
	int next=(rank+1)%size;
	int prev=(rank+size-1)%size;

	int i;
	if (rank%2==0)
	  {
		send(context().peer(next),rank); 
		std::cout << rank << std::endl;
		recv(context().peer(prev),i);
	  }
	else 
	  {
		recv(context().peer(prev),i);
		std::cout << rank << std::endl;
		send(context().peer(next),rank); 
	  }
  }
};


int main (int argc, char *argv[])
{
  TT::ThreadPool<MyThread> pool(10,MyThread());
  pool.start();
  pool.stop();
}
