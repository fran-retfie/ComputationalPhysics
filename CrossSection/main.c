
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"
#include <RecFormulae.h>
#define N 10000
#define Lmax 5
#define Nmax 3

const double tEnd = 2;
const double h = tEnd/(double) N;
const double Estep = 0.01;
double sigma = 3.18e-10; //angstrom
double epsilon = 5.9;//meV
const double h2m = 0.03517; //hbar/2m in units of sigma and epsilon = 7.77x10^(-19)
const double b10 = 8./(25.*h2m);
//double b10 = sqrt(arg); //b^5 in units of sigma and epsilon
const double xStart = 0.7;

double delta_l[Lmax];
double Ea[2][Nmax*Lmax];
double psi[2][N];
double psi_norm[2][N];
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
  return (-4/h2m*(x12 - x6) - l*(l+1)*x2 + Ex);//
}
double initialCondition(double r0)
{
  double x = 1./r0;
  double x2 = x*x;
  double x5 = x2*x2*x;
  return exp(-sqrt(b10)*x5);
}
double integrate(double Ex,int l)
{
  double (*f)(double,double) = &fr;

  //initial condition
  psi[0][0] = xStart; //in units of sigma r_low is in [0,1]
  a[0] = initialCondition(psi[0][0]);
  psi[1][0] = a[0];

  psi[0][1] = xStart+h; //in units of sigma r_low is in [0,1]
  a[1] = initialCondition(psi[0][1]);
  psi[1][1] = a[1];

  psi[0][2] = xStart+h*2; //in units of sigma r_low is in [0,1]
  a[2] = initialCondition(psi[0][2]);
  psi[1][2] = a[2];

  for(int j = 3; j < N; j++)
  {
    psi[0][j] = h*j + xStart;
    a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];
  }
  return (a[N-1]-a[N-3]);//
  // //initial condition
  // psi[0][0] = 0.8*sigma;
  // psi[0][1] = psi[0][0]+h; //in units of sigma r_low is in [0,1]
  // double x = 1./psi[0][0];
  // double x2 = x*x;
  // double x5 = x2*x2*x;
  // a[0] = exp(-sqrt(b10)*x5);
  // a[1] = pow(a[0],(double)l);
  //
  // for(int j = 2; j < N; j++)
  // {
  //   psi[0][j] = psi[0][0] + h*j;
  //   a[j] = NumerovInt(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
  //   psi[1][j] = a[j]/psi[0][j];
  // }
  //
  // psi[1][0] = psi[1][2];
  // psi[1][1] = psi[1][2];

  // return a[N-1];//
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
  epsilon = 1;
  //parameters to determine change in sign of solution
  int ee=0;
  for(E = epsilon*0.01; E<3.5; E += Estep)
  {
    ee++;
  }
  double sigmaTOT[ee];
  double k,x;
  double K;
  double *j1,*j2;
  double *n1,*n2;
  double s[Lmax];
  double r1,r2;
  double R1[Lmax],R2[Lmax];
  double tmp0;
  char filename[17];
  r1 = 20;//in units of sigma
  r2 = 50;
  int in=0;
  double sum=0;
  double pi = 4.*atan(1.);
  ee=0;
  for(E = epsilon*0.01; E<3.5; E += Estep)
  {

    //scan on angular momentum
    for(l = 0; l<Lmax; l++)
    {
      tmp0 = integrate(E,l);
      double norm = CalcNorm();
      for(int j = 0; j < N; j++)
      {
        psi_norm[1][j] = psi[1][j]/norm;
        psi_norm[0][j] = psi[0][j];
        if ((psi[0][j] > 20)&&(psi[0][j]<21))
        {
          r1 = psi[0][j];
          R1[l] = psi[1][j];
        }
        if ((psi[0][j] > 50)&&(psi[0][j]<51))
        {
          r2 = psi[0][j];
          R2[l] = psi[1][j];
        }
      }
      // printf("r1=%g    r2=%g\n", r1,r2);
      sprintf(filename,"dati/plot%02i_%02i",l,ee);
      writeCSVdouble(filename,(double *)psi_norm, 2, N);
    }
    ee++;
    k = sqrt(E/h2m);
    //calculate bessel function for given E
    x = k*r1;
    // j1 = Bessel(x,Lmax,0,s);
    // n1 = Bessel(x,Lmax,1,s);
    // x = k*r2;
    // j2 = Bessel(x,Lmax,0,s);
    // n2 = Bessel(x,Lmax,1,s);
    // //calculate delta_l
    // for(l = 0; l<Lmax; l++)
    // {
    //   K = (R1[l]*r2)/(R2[l]*r1);
    //   delta_l[l] = atan((K*j2[l] - j1[l])/(K*n2[l] - n1[l]));
    //   sum = sum + (2*l+1)*sin(delta_l[l])*sin(delta_l[l]);
    // }
    // sigmaTOT[in] = 4.*pi/(k*k)*sum*sigma*sigma;
    //printf("E = %g,    sigmaTOT_%i = %g\n",E, in, sigmaTOT[in]);
    sum = 0;
    in++;
  }


  writeCSVdouble("energies.csv",(double *)Ea, 2, Nmax*Lmax);
  return 0;
}
