
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000
#define Nmax 5

const double tEnd = 10;
const double h = tEnd/(double) N;
const double Estep = 0.001;
double Ea[2][Nmax];
double psi[2][N];
double a[N];

double E = 0;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  double x2 = x*x;
  return (-x2*0.25 + Ex);
}

double integrate(double Ex, int n)
{
  double (*f)(double,double) = &fr;

  //initial condition
  if(n%2==1)
  {
    psi[0][0] = 0;
    a[0] = 0;
    psi[0][1] = h;
    a[1] = h;
  }
  else
  {
    psi[0][0] = 0;
    a[0] = 1;
    psi[0][1] = h;
    a[1] = 1;
  }

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

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

    E = Estep;
    //scan on energy
    for(int n = 0; n<Nmax; n++)
    {
      Eold = E;

    //simple scan
    for(; (gamma*gammaOld)>0|(E < Eold+Estep*3); E += Estep)
    {
      gammaOld = gamma;
      gamma = integrate(E,n);
    }

    //tangent method
    double E1 = E - Estep;
    double E2;
    for(int j = 0; j<100; j++)
    {
      double g = integrate(E,n);
      double g1 = integrate(E1,n);
      if(g==g1) break;
      E2 = E;
      E = (E1*g-E*g1)/(g - g1);
      if(fabs(E-E1)<1e-14*Estep) break;
      E1 = E2;
    }

    printf("n = %i, E = %15.14g\n",n,E);
    Ea[1][n] = E;

    double norm = CalcNorm();
    for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

    char filename[11];
    sprintf(filename,"dati/plot%02i",n);
    writeCSVdouble(filename,(double *)psi, 2, N);

  }

  writeCSVdouble("energies.csv",(double *)Ea, 2, Nmax);
  return 0;
}
