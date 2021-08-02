#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

double alpha[4] = {14.899983, 2.726485, 0.757447, 0.251390};

double Chi(double r, double a)
{
  double r2 = (r-1)*(r-1);//1 va sostituito con distanza media orbitale

  return (exp(-a*r2));
}

double Vc(double r1, double r2)
{
  return 1./fabs(r1-r2);
}

double f_exch(double r1, double r2, int p, int q, int r, int s)
{
  return Chi(r1,alpha[p])*Chi(r2,alpha[r])*Vc(r1,r2)*Chi(r1,alpha[s])*Chi(r2,alpha[q]);
}

double doubleInt(double f)
{

}

int main()
{
  //Compute function at different positions storing each value
  //---check where function goes to zero (below threshold)
  double threshold = 1e-5;
  double check = 1;
  int n1,n2;
  int p,q,r,s;
  double r1,r2,h;
  r1 = 1;
  r2 = 1;
  h = 1e-3;



  for(p=0;p<4;p++)
  {
    for(q=0;q<4;q++)
    {
      for(r=0;r<4;r++)
      {
        for(s=0;s<4;s++)
        {
          check = 1;
          n1 = 0;
          r1 = 1;
          r2 = 1;
          //loop on r1
          while(check>threshold)
          {
            check = f_exch(r1,r2,p,q,r,s);
            r1 += h;
            n1++;
          }
          printf("-------------[p,q,r,s] = [%d,%d,%d,%d]-------------\n",p,q,r,s );
          printf("n = %d      r1 = %e\n", n1, r1);
          check = 1;
          n2 = 0;
          r1 = 1;
          r2 = 1;
          while(check>threshold)
          {
            check = f_exch(r1,r2,p,q,r,s);
            r2 += h;
            n2++;
          }
          printf("n = %d      r2 = %e\n", n2, r2);
        }
      }
    }
  }



  return 0;
}
