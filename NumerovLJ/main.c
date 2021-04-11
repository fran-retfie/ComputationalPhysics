
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 1000
#define Lmax 3
#define Nmax 3

const double Estep = 0.001;
//const double sigma = 3.18; //angstrom
const double epsilon = 5.9; //<-normalized-energy----da-aggiustare-il-valore-----------
const double h2m = 1; //hbar/2m in units of sigma and epsilon = 7.77x10^(-19)
const double b10 = 8./(25.*h2m);
//double b10 = sqrt(arg); //b^5 in units of sigma and epsilon
const double xStart = 0.7;
const double tEnd = 4;
const double h = tEnd/(double) N;

double Ea[2][Nmax*Lmax];
double psi[2][N];
double a[N];

int l = 0;
double E;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  double x2 = 1./(x*x);
  double x4 = x2*x2;
  double x6 = x4*x2;
  double x12 = x6*x6;
  return (-4*epsilon*(x12 - x6) - l*(l+1)*x2 + Ex);//
}

double initialCondition(double r0)
{
  double x = 1./r0;
  double x2 = x*x;
  double x5 = x2*x2*x;
  return exp(-sqrt(b10)*x5);
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;

  //initial condition
  psi[0][0] = xStart; //in units of sigma r_low is in [0,1]
  a[0] = initialCondition(psi[0][0]);
  psi[1][0] = a[0]/psi[0][0];

  psi[0][1] = xStart+h; //in units of sigma r_low is in [0,1]
  a[1] = initialCondition(psi[0][1]);
  psi[1][1] = a[1]/psi[0][1];

  psi[0][2] = xStart+h*2; //in units of sigma r_low is in [0,1]
  a[2] = initialCondition(psi[0][2]);
  psi[1][2] = a[2]/psi[0][2];

  for(int j = 3; j < N; j++)
  {
    psi[0][j] = h*j + xStart;
    a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j]/psi[0][j];
  }

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
  printf("%f\n", b10);
  psi[1][N-1] = 1;

  double gamma = 1;
  double gammaOld = 1;
  double Eold;

  //scan on angular momentum
  for(l = 0; l<Lmax; l++)
  {
    E = -epsilon;

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
      for(int j = 0; j<100; j++)
      {
        double g = integrate(E);
        double g1 = integrate(E1);
        if(g==g1) break;
        E2 = E;
        E = (E1*g-E*g1)/(g - g1);
        if(fabs(E - E1)<1e-17) break;
        E1 = E2;
      }

      printf("l = %i, n = %i, E = %15.14g\n",l,n,E);
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
