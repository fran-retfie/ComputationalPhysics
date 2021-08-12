#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include "CSVio.h"

double alpha[4] = {14.899983, 2.726485, 0.757447, 0.251390};
int n = 4;
double Chi(double r, double a)
{
  double r2 = r*r;//1 va sostituito con distanza media orbitale

  return (exp(-a*r2));
}

double Vc(double r1, double r2)
{
  return 1./fabs(r1-r2);
}

double f_exch(double r1, double r2, int p, int q, int r, int s)
{
  return Chi(r1,alpha[p])*Chi(r2,alpha[r])*Vc(r1,r2)*Chi(r1,alpha[s])*Chi(r2,alpha[q]);
  // return r1*r1+r2*r2;
}

double E_pq(int p, int q, double h)
{
   double threshold = 1e-5;
   double check = 1;

   double r1,r2,integ_rs,integ_r,integ;
   int s,r,i,j,n1,n2;
   for(r=0;r<n;r++)
   {
     for(s=0;s<n;s++)
     {
       check = 1;
       n1 = 0;
       r1 = 1;
       r2 = 1;
       while(check>threshold)
       {
         check = f_exch(r1,r2,p,q,r,s);
         r1 += h;
         n1++;
       }
       n1 += (1/h);
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
       n2 += (1/h);
     }
   }
   double fx[n1];
   for(r=0;r<n;r++)
   {
     for(s=0;s<n;s++)
     {
       for(i=0;i<n1;i++)
       {
         fx[i] = 0;
         r1 = h + h*i;
         for(j=0;j<n2;j++)
         {
           r2 = h + h*j;
           if(r1==r2)
           {
             fx[i] += 0;
           }
           else
           {
             if(j==0||j==n2-1) fx[i] += f_exch(r1,r2,p,q,r,s);
             else if(j%2==0)   fx[i] += 2.*f_exch(r1,r2,p,q,r,s);
             else              fx[i] += 4.*f_exch(r1,r2,p,q,r,s);
             // if(r1==r2) printf("%d %d %f r = %f\n", i,j, fx[i],r2);
             // printf("%d\n", j);
           }
         }
         fx[i] *= (h/3.);
         // printf("%d   %f\n", i, fx[i]);
       }
       integ_rs = 0;
       for(i=0;i<n1;i++)
       {
         if(i==0||i==n1-1) integ_rs += fx[i];
         else if(i%2==0)   integ_rs += 2.*fx[i];
         else              integ_rs += 4.*fx[i];
         // printf("%d\n", i);
       }
       integ_rs *= (h/3);
       // printf("I_rs=%f, s=%d,r=%d\n", integ_rs,s,r);
       integ += integ_rs;
     }

   }
   printf("[%d,%d]-------integral = %f-------\n", p,q,integ);
   return integ;
}

double E_pqrs(int p, int q, int r, int s, double h)
{
   double threshold = 1e-5;
   double check = 1;

   double r1,r2,integ_rs,integ_r,integ;
   int i,j,n1,n2;

   check = 1;n1 = 0;r1 = 1;r2 = 1;
   while(check>threshold)
   {
     check = f_exch(r1,r2,p,q,r,s);
     r1 += h;
     n1++;
   }
   n1 += (1/h);
   check = 1;n2 = 0;r1 = 1;r2 = 1;
   while(check>threshold)
   {
     check = f_exch(r1,r2,p,q,r,s);
     r2 += h;
     n2++;
   }
   n2 += (1/h);

   double fx[n1];

   for(i=0;i<n1;i++)
   {
     fx[i] = 0;
     r1 = h + h*i;
     for(j=0;j<n2;j++)
     {
       r2 = h + h*j;
       if(r1==r2)
       {
         fx[i] += 0;
       }
       else
       {
         if(j==0||j==n2-1) fx[i] += f_exch(r1,r2,p,q,r,s);
         else if(j%2==0)   fx[i] += 2.*f_exch(r1,r2,p,q,r,s);
         else              fx[i] += 4.*f_exch(r1,r2,p,q,r,s);
         // if(r1==r2) printf("%d %d %f r = %f\n", i,j, fx[i],r2);
         // printf("%d\n", j);
       }
     }
     fx[i] *= (h/3.);
     // printf("%d   %f\n", i, fx[i]);
   }
   integ_rs = 0;
   for(i=0;i<n1;i++)
   {
     if(i==0||i==n1-1) integ_rs += fx[i];
     else if(i%2==0)   integ_rs += 2.*fx[i];
     else              integ_rs += 4.*fx[i];
     // printf("%d\n", i);
   }
   integ_rs *= (h/3);
   // printf("-----I_rs=%f -- s=%d,r=%d------\n", integ_rs,s,r);

   return integ_rs;
}


int main()
{
  //Compute function at different positions storing each value
  //---check where function goes to zero (below threshold)
  double threshold = 1e-5;
  double check = 1;
  int n1,n2;
  int p,q,r,s;
  double h;
  double Exch_pqrs[n][n];
  double Dir_pqrs[n][n];
  h = 1e-3;

  char filename[50];

  for(p=0;p<n;p++)
  {
    printf("---- p = %d ----\n", p);
    for(q=0;q<n;q++)
    {
      for(r=0;r<n;r++)
      {
        for(s=0;s<n;s++)
        {
          Exch_pqrs[r][s] = -E_pqrs(p,q,r,s,h);
          Dir_pqrs[r][s] = E_pqrs(p,s,r,q,h);
        }
      }
      sprintf(filename,"Exc/p%02iq%02i.csv", p,q);
      writeCSVdouble(filename, (double *)Exch_pqrs, n, n);
      sprintf(filename,"Dir/p%02iq%02i.csv", p,q);
      writeCSVdouble(filename, (double *)Dir_pqrs, n, n);
    }
  }


  // OLD: calcola E_pq (a me serve E_pqrs)
  // double Exch_pq[n][n];
  //
  // for(p=0;p<4;p++)
  // {
  //   for(q=0;q<4;q++)
  //   {
  //     Exch_pq[p][q] =  E_pq(p,q,h);
  //   }
  // }
  //
  // writeCSVdouble("Ex_terms.csv", (double *)Exch_pq, n, n);



  return 0;
}
