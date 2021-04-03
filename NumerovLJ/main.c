
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000
#define Lmax 3
#define Nmax 3

const double tEnd = 500;
const double h = tEnd/(double) N;
const double Estep = 0.001;
const double sigma = 3.18; //angstrom
const double epsilon = 5.9;//meV
const double h2m = 1; //hbar/2m in units of sigma and epsilon = 7.77x10^(-19)
const double b10 = 4./(25.*h2m);
//double b10 = sqrt(arg); //b^5 in units of sigma and epsilon



double Ea[2][Nmax*Lmax];
double psi[2][N];
double a[N];

int l = 0;
double E = 0;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  double x2 = 1./(x*x);
  double x4 = x2*x2;
  double x6 = x4*x2;
  double x12 = x6*x6;
  return (-4*(x12 - x6) - l*(l+1)*x2 + Ex);//
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;

  //initial condition
  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = 0.5; //in units of sigma r_low is in [0,1]
  double x = 1./psi[0][1];
  double x2 = x*x;
  double x5 = x2*x2*x;
  a[1] = exp(-sqrt(b10)*x5);

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt(psi[0][j], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j]/psi[0][j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

  return (a[N-1]-a[N-3]);//
}

int main()
{
  printf("%f\n", b10);
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

      printf("l = %i, n = %i, E = %15.14g\n",l,n,E);
      Ea[1][n+l*Nmax] = E;
      Ea[0][n+l*Nmax] = (double) l;

      char filename[13];
      sprintf(filename,"dati/plot%02i%02i",n,l);
      writeCSVdouble(filename,(double *)psi, 2, N);
    }
  }


  writeCSVdouble("energies.csv",(double *)Ea, 2, Nmax*Lmax);
  return 0;
}
