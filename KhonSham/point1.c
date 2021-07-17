
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000

//max number of shells
#define Nmax 4
#define double_t __float128


const double_t tEnd = 30;
const double_t h = tEnd/(double_t) N;
const double_t Estep = 0.001;
double_t psi[3][N];
double_t a[N];
double_t Ea[2][Nmax];

double_t E = 0;

//Na atom
const double_t rs = 3.93;
//K atom
//const double_t rs = 4.86;
const double_t rs3 = (rs*rs*rs);
const double_t o_rs3 = 1/rs3;
double_t Rc;




//effective potential (normalized!!)
double_t fr(double_t r, double_t Ex)
{
  double_t Vext = o_rs3 * ( (r < Rc) ? (r*r - 3*Rc*Rc) : (-2*Rc*Rc*Rc/r) );
  return (-Vext + Ex);
}

double_t integrate(double_t Ex)
{
  double_t (*f)(double_t,double_t) = &fr;

  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = h;
  a[1] = h;
  psi[0][2] = 2*h;
  a[2] = 2*h;

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt_t(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

  return (a[N-1]-a[N-2]);
}

double_t CalcNorm()
{
  double_t I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += h/3*(psi[1][j-1]*psi[1][j-1]+4*psi[1][j]*psi[1][j]+psi[1][j+1]*psi[1][j+1]);
  }
  return sqrtq(I);
}

int main()
{
  for(int n = 1; n <= Nmax; n++)
  {
    int Ne = 0;
    for(int i=0; i<n; i++) Ne += 2*(2*i+1);

    Rc = cbrtq(Ne)*rs;

    E = -10;

    double_t gamma = 1;
    double_t gammaOld = 1;
    double_t Eold = E;

    for(; (gamma*gammaOld)>0|(E < Eold+Estep*4); E += Estep)
    {
      gammaOld = gamma;
      gamma = integrate(E);
      if(E>0) break;
    }

    double_t E1 = E - 2*Estep;
    double_t E2;
    double_t g;
    double_t g1;

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
      if(fabsq(g-g1)<1e-100) break;
      if(j>200) printf("ooo\n");
    }

    for(int j = 0; j < N; j++) psi[2][j] = fr(j*h,0);

    printf("N_e = %i,\t E = %Qe\n",Ne,E);
    Ea[1][n] = E1;

    double_t norm = CalcNorm();
    for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

    char filename[11];
    sprintf(filename,"dati/plot%02i",n);
    writeCSVdouble_t(filename,(double_t *)psi, 3, N);
  }

  writeCSVdouble_t("dati/energies.csv",(double_t *)Ea, 2, Nmax);
  return 0;
}
