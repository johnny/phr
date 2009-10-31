#ifndef TIMER_H
#define TIMER_H

// headers for getrusage(2)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

// headers for stderror(3)
#include <string.h>

// access to errno
#include <errno.h>

// printf etc.
#include <stdio.h>

/* a simple stop watch

This functions report the elapsed user-time, i.e. time spent computing,
after the last call to Timer::reset(). The results are seconds and
fractional seconds. Note that the resolution of the timing depends
on your OS kernel which should be somewhere in the milisecond range.

The functions are basically wrappers for the libc-function getrusage().

*/

//! reset timer
inline
void reset_timer(struct timeval * timer)
{
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru))
    fprintf(stderr, "Error: %s\n", strerror(errno));
  *timer = ru.ru_utime;
}

//! get elapsed user-time in seconds
inline
double get_timer(struct timeval timer)
{
  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru))
    fprintf(stderr, "Error: %s\n", strerror(errno));
  return 1.0 * (ru.ru_utime.tv_sec - timer.tv_sec) + (ru.ru_utime.tv_usec - timer.tv_usec) / (1000.0 * 1000.0);
}

//! print elapsed user-time in seconds
inline
void print_timer(struct timeval timer)
{
  printf("elapsed time: %f s\n", get_timer(timer));
}

#endif
