#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>

double get_time ()
{
/*   struct rusage ru; */
/*   struct timeval cstop; */
  
/*   getrusage(RUSAGE_SELF, &ru); */
/*   cstop = ru.ru_utime; */
/*   return 1.0*cstop.tv_sec + cstop.tv_usec/1000000.0; */
/*   return clock()/CLOCKS_PER_SEC; */
  return omp_get_wtime();
}
