
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 400

const double tEnd = 50;
const double h = tEnd/(double) N;
const double Estep = 0.01;
double psi[2][N];
double a[N];

int l = 0;
double E = 0;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  return (-(x*x)/9/4+E);
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;

  //initial condition
  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = h;
  a[1] = h;

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt(psi[0][j], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j]/psi[0][j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

  return (psi[1][N-1]-psi[1][N-3]);//psi[1][N-2];
}

int main()
{
  psi[1][N-1] = 1;

  double gamma = 1;
  double gammaOld = 1;
  double Eold;

  for(int i = 0; i<10; i++)
  {
    Eold = E;
    for(; (gamma*gammaOld)>0|(E < Eold+Estep*3); E += Estep)
    {
      gammaOld = gamma;
      gamma = integrate(E);
    }

    double E1 = E - Estep;
    double E2;
    for(int j = 0; j<1; j++)
    {
      E2 = E;
      E = (E1*integrate(E)-E*integrate(E1))/(integrate(E) - integrate(E1));
      E1 = E2;
    }

  printf("%lg\n", E);
  }

  writeCSVdouble("dati.csv",(double *)psi, 2, N);
  return 0;
}
