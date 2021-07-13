
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000

//max number of shells
#define Nmax 4

const double tEnd = 23;
const double h = tEnd/(double) N;
const double Estep = 0.001;
double psi[3][N];
double a[N];
double Ea[2][Nmax];

double E = 0;

//Na atom
const double rs = 3.93;
//K atom
//const double rs = 4.86;
const double rs3 = (rs*rs*rs);
const double o_rs3 = 1/rs3;
double Rc;




//effective potential (normalized!!)
double fr(double r, double Ex)
{
  double Vext = o_rs3 * ( (r < Rc) ? (r*r - 3*Rc*Rc) : (-2*Rc*Rc*Rc/r) );
  return (-Vext + Ex);
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;

  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = h;
  a[1] = h;
  psi[0][2] = 2*h;
  a[2] = 2*h;

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

  return (a[N-1]-a[N-2]);
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
  for(int n = 1; n <= Nmax; n++)
  {
    int Ne = 0;
    for(int i=0; i<n; i++) Ne += 2*(2*i+1);

    Rc = cbrt(Ne)*rs;

    E = -10;

    double gamma = 1;
    double gammaOld = 1;
    double Eold = E;

    //simple scan
    for(; (gamma*gammaOld)>0|(E < Eold+Estep*3); E += Estep)
    {
      gammaOld = gamma;
      gamma = integrate(E);
    }

    //tangent method
    double E1 = E - Estep;
    double E2;
    for(int j = 0; j<100; j++)
    {
      double g = integrate(E);
      double g1 = integrate(E1);
      if(g==g1) break;
      E2 = E;
      E = (E1*g-E*g1)/(g - g1);
      if(fabs(E-E1)<1e-15*Estep) break;
      E1 = E2;
    }

    for(int j = 0; j < N; j++) psi[2][j] = fr(j*h,0);

    printf("N_e = %i,\t E = %15.14g\n",Ne,E);
    Ea[1][n] = E1;

    double norm = CalcNorm();
    for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

    char filename[11];
    sprintf(filename,"dati/plot%02i",n);
    writeCSVdouble(filename,(double *)psi, 3, N);
  }

  writeCSVdouble("dati/energies.csv",(double *)Ea, 2, Nmax);
  return 0;
}
