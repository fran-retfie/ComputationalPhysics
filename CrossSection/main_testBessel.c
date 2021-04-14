#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <RecFormulae.h>

int main()
{
  double E,r,h2m;
  int Lmax,kind;
  E = 1.5;
  r = 5.4;
  h2m = 1;
  Lmax = 10;
  kind = 1;
  double s[Lmax];
  double *k;
  k = Bessel(E,r,Lmax,h2m,kind,s);
  for(int l=0;l<Lmax;l++)
  {
    printf("k[%i] = %g\n",l,k[l]);
  }
  return 0;
}
