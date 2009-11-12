#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

// check primeness of n, brute force method
inline bool is_prime(int n)
{
  const int i = (int) std::sqrt(n);
  for (int k=2; k<i+1; ++k)
    if (n%k == 0)
      return false;

  return true;
}

// check Goldbach conjecture by counting number pairs (i,j)
// satisfying it (i+j == n)
int goldbach_pairs(int n)
{
  if ((n+1)%2 == 0) return 0;   // odd numbers not covered by Goldbach
  int result = 0;               // number of number pairs (i,j) with (i+j) == n

  for (int i=2; i<(n/2+1); ++i) // loop through all i = 2...n/2
    if (is_prime(i))            // check primeness of i
      if (is_prime(n-i))        // check primeness of of n-i
          ++result;

  return result;
}

//-------------------------------------------------------------------------
// main -- computes number of goldbach pairs between 1 and N
//-------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // check program arguments
  if (argc != 2)
  {
    std::cout << "Usage: ./<executable> N" << std::endl;
    exit (EXIT_FAILURE);
  }
  int N = atoi(argv[1]);

  // vector storing number of found pairs for each number
  std::vector<int> goldbach(N, 0);

  // find and count Goldbach pairs
  for (int i=1; i<=N; i++)
    goldbach[i-1] = goldbach_pairs(i);

  // output number of Goldbach pairs
  for (int i=0; i<goldbach.size(); ++i)
    std::cout << i+1 << " Number of Goldbach pairs:  " << goldbach[i] << std::endl;

  return 0;
}
