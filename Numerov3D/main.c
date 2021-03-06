
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000
#define Lmax 3
#define Nmax 3

const double tEnd = 10;
const double h = tEnd/(double) N;
const double Estep = 0.001;
double Ea[2][Nmax*Lmax];
double psi[2][N];
double a[N];

int l = 0;
double E = 0;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  if(x == 0) printf("error\n");
  double x2 = x*x;
  return (-x2*0.25 - l*(l+1)/x2 + Ex);
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;

  //initial condition
  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = h;
  a[1] = h;
  psi[0][2] = 2*h;
  a[2] = 2*h;

  for(int j = 3; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j]/psi[0][j];
  }

  psi[1][2] = a[2]/psi[0][2];
  psi[1][1] = a[1]/psi[0][1];
  psi[1][0] = psi[1][1];

  return (a[N-1]-a[N-3]);//
}

double CalcNorm()
{
  double I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += h/3*(psi[1][j-1]*psi[1][j-1]+4*psi[1][j]*psi[1][j]+psi[1][j+1]*psi[1][j+1]);
  }
  return sqrt(I);
}

int main()
{
  psi[1][N-1] = 1;

  double gamma = 1;
  double gammaOld = 1;
  double Eold;

  //scan on angular momentum
  for(l = 0; l<Lmax; l++)
  {
    E = Estep;
    //scan on energy
    for(int n = 0; n<Nmax; n++)
    {
      Eold = E;

      //simple scan
      for(; (gamma*gammaOld)>0|(E < Eold+Estep*3); E += Estep)
      {
        gammaOld = gamma;
        gamma = integrate(E);
      }

      //tangent method
      double E1 = E - Estep;
      double E2;
      for(int j = 0; j<8; j++)
      {
        double g = integrate(E);
        double g1 = integrate(E1);
        if(g==g1) break;
        E2 = E;
        E = (E1*g-E*g1)/(g - g1);
        E1 = E2;
      }

      printf(" %i & %i & %15.14g & %3.1g \\\\ \n",l,n,E,(E-(2*n+l+1.5))/(2*n+l+1.5));
      Ea[1][n+l*Nmax] = E;
      Ea[0][n+l*Nmax] = (double) l;

      double norm = CalcNorm();
      for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

      char filename[13];
      sprintf(filename,"dati/plot%02i%02i",n,l);
      writeCSVdouble(filename,(double *)psi, 2, N);
    }
  }

  writeCSVdouble("dati/energies.csv",(double *)Ea, 2, Nmax*Lmax);
  return 0;
}
