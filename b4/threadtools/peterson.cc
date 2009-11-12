#include<iostream>
#include"tt.hh"

int in[2]={0,0}, last=0;

class Process : public TT::BasicThread {
public:
  Process (int i) : rank(i) {}
  virtual void run () {
	in[rank] = 1;
	last = rank;
	while (in[1-rank]==1 && last==rank) ;
	std::cout << "Process " << rank << " in critical section" << std::endl;
	in[rank] = 0;
  }
private:
  int rank;
};

int main (int argc, char *argv[])
{
  Process p0(0),p1(1);
  p0.start(); p1.start();
  p0.stop(); p1.stop();
}
