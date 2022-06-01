#include "stdio.h"
#include "RecFormulae.h"
#include "math.h"
#include "stdlib.h"
void Bessel(double x, int Lmax, int kind, double s[])
{
//  x=300;
  if(kind==0) //j_l
  {
    s[1] = sin(x)/x; //s_(0)
    s[0] = cos(x)/x; //s_(-1)
    for(int l=0; l<Lmax;l++)
    {
      s[l+2] = (2*l+1)/x*s[l+1] - s[l];
    }
  }
  else //n_l
  {
    s[1] = -cos(x)/x; //s_(0)
    s[0] = sin(x)/x; //s_(-1)
    for(int l=0; l<Lmax;l++)
    {
      s[l+2] = (2*l+1)/x*s[l+1] - s[l];
    }
  }
  for(int l=0; l<Lmax+1;l++) s[l] = s[l+1];
}
