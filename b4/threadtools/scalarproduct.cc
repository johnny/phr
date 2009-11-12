#include"tt.hh"

// shared memory scalar product 
const int P=10;
const int N=100000;
double x[N],y[N];
double s=0;

class Process : public TT::ThreadPoolElement<Process>
{
public:
  virtual void run ()
  {
	const int p = context().rank();
	for (int i=p*N/P; i<(p+1)*N/P; i++)
	  s += x[i]*y[i];   // Vorsicht!
  }
};

int main (int argc, char *argv[])
{
  TT::ThreadPool<Process> pool(10,Process());
  pool.start(); pool.stop();
}
