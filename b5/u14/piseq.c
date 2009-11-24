#include <stdio.h>
#include <math.h>

/* Berechnung von pi ueber den arcus-tangens */
int main(int argc, char **argv)
{
  const double PI25DT = 3.141592653589793238462643;
  long n, i;
  double pi, h, sum, x;

  while (n!=0)
  {
    printf("Anzahl Intervalle (0 beendet): ");
    scanf("%d",&n);

    if (n == 0)
      break;
    else
    {
      h = 1.0/(double) n;
      sum = 0.0;
      for (i=1; i<=n; i+=1)
      {
        x = h * ((double) i - 0.5);
        sum += (4.0 / (1.0 + x*x));
      }
      pi = h * sum;

      printf("pi: %.16f, Fehler: %.16f\n",
          pi, fabs(pi - PI25DT));
    }
  }

  return 0;
}
