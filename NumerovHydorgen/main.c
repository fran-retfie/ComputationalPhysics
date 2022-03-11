
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"
#include <gsl/gsl_sf_exp.h>

#ifndef long_double_t
  typedef  long double long_double_t; //__float128  long_double_t;
#endif

#define N 1000
#define Lmax 4
#define Nmax 5

const long_double_t tStart = 0.01;
const long_double_t tEnd = 10;
const long_double_t h = (tEnd-tStart)/(long_double_t) N;
const long_double_t Estep = 0.01;
long_double_t Ea[2][Nmax*Lmax];
long_double_t psi[3][N];
long_double_t a[N];

int n = 1;
int l = 0;
long_double_t E = 0;

const double epsilon = 10;
const double v0 = 27.4080280210014; // epsilon = 1
const double aaa = 7.39;
const double bbb = 3.22;
const double u1 = 8.95644589285677; // epsilon = 1
const double u2 = 4.478222946428356; // epsilon = 1

const long_double_t C = epsilon/4; //13.605693122904;

//effective potential (normalized!!)
long_double_t fr(long_double_t x, long_double_t Ex)
{
  //if(x == 0) printf("error\n");
  //long_double_t r1 = 1/x;
  //long_double_t pot = (2*r1 + Ex);
  //if(l != 0) pot -= l*(l+1)*r1*r1;
  double pot = + v0/x * (aaa*gsl_sf_exp(-u1*x) - bbb*gsl_sf_exp(-u2*x)) + Ex;
  return pot;
}

long_double_t integrate(long_double_t Ex)
{
  long_double_t (*f)(long_double_t,long_double_t) = &fr;

  //initial condition
  psi[0][0] = tStart;
  a[0] = -pow(-psi[0][0],l+1)*exp(-psi[0][0]/n);
  psi[0][1] = h + tStart;
  a[1] = -pow(-psi[0][1],l+1)*exp(-psi[0][1]/n);
  psi[0][2] = 2*h + tStart;
  a[2] = -pow(-psi[0][2],l+1)*exp(-psi[0][2]/n);

  for(int j = 3; j < N; j++)
  {
    psi[0][j] = h*j  + tStart;
    a[j] = NumerovIntLD(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];//psi[0][j];
    psi[2][j] = -fr(psi[0][j],E);
  }

  psi[1][2] = a[2];//psi[0][2];
  psi[1][1] = a[1];//psi[0][1];
  psi[1][0] = psi[1][1];

  return (a[N-1]-a[N-3]); //a[N-1]-exp(-psi[0][N-1]/n);
}

long_double_t CalcNorm()
{
  long_double_t I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += h/3*(psi[1][j-1]*psi[1][j-1]+4*psi[1][j]*psi[1][j]+psi[1][j+1]*psi[1][j+1]);
  }
  return sqrt(I);
}

int main()
{
  printf("%i\n", sizeof(long_double_t));

  psi[1][N-1] = 1;

  long_double_t gamma = 1;
  long_double_t gammaOld = 1;
  long_double_t Eold;

  //scan on angular momentum
  for(l = 0; l<Lmax; l++)
  {
    E =-2;

    //scan on energy
    for(n = 1; n<Nmax; n++)
    {
      Eold = E;
      //simple scan
      for(; (gamma*gammaOld)>0|(E < Eold+Estep*4); E += Estep)
      {
        gammaOld = gamma;
        gamma = integrate(E);
        if(E>0) break;
      }


      long_double_t E1 = E - 2*Estep;
      long_double_t E2;
      long_double_t g;
      long_double_t g1;

/*
      //tangent method
      for(int j = 0; j<500; j++)
      {
        g = integrate(E);
        g1 = integrate(E1);
        //if(j>200) printf("%g\n",g);
        if(g==g1) break;
        E2 = E;
        E = (E1*g-E*g1)/(g - g1);
        E1 = E2;
        //if(fabsl(g-g1)<1e-200) break;
      }
*/

      //bisection
      for(int j = 0; j<200; j++)
      {
        E2 = (E+E1)*0.5;
        g = integrate(E);
        g1 = integrate(E2);
        if(g*g1<0)
        {
          //printf("O");
          E1 = E2;
        }
        else
        {
          //printf("X");
          E = E2;
        }
        //if(fabs(g-g1)<1e-100) break;
        //if(j>200) printf("%g\n",E);
      }
      printf("\n");


      printf("l = %i, n = %i, E = %40.39Lg, Err = %40.39Lg\n",l,n, (E*C), (E+1));
      Ea[1][n+l*Nmax] = E*C;
      Ea[0][n+l*Nmax] = (long_double_t) l;

      long_double_t norm = CalcNorm();
      for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

      char filename[13];
      sprintf(filename,"dati/plot%02i%02i",n,l);
      writeCSVlongdouble(filename,(long_double_t *)psi, 3, N);

    }
  }

  writeCSVlongdouble("energies.csv",(long_double_t *)Ea, 2, Nmax*Lmax);
  return 0;
}
