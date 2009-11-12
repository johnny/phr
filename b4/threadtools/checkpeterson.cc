#include<iostream>
#include"tt.hh"

int in[2]={0,0}, last=0;
int count=0;
const int sections=100000000;

class Process : public TT::BasicThread {
public:
  Process (int i) : rank(i) {}
  virtual void run () {
	for (int i=0; i<sections; i++) {
 	  in[rank] = 1;
 	  last = rank;
 	  while (in[1-rank]==1 && last==rank) ;
	  count += 1;
 	  in[rank] = 0;
	}
  }
private:
  int rank;
};

int main (int argc, char *argv[])
{
  Process p0(0),p1(1);
  p0.start(); p1.start();
  p0.stop(); p1.stop();
  std::cout << "    count = " << count << std::endl;
  std::cout	<< " sections = " << 2*sections << std::endl;
  std::cout	<< "   failed = " << 2*sections-count << std::endl;
}
