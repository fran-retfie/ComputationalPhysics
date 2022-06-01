
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"
#include <RecFormulae.h>
#define N 10000
#define Lmax 7
#define Nmax 3
#define sizeEa 551
double Emax;
const double tEnd = 100;
const double h = tEnd/(double) N;
const double Estep = 0.001;
double sigma = 3.18e-10; //angstrom
double epsilon = 5.9;//meV
const double h2m = 0.03517; //hbar/2m in units of sigma and epsilon = 7.77x10^(-19)
const double b10 = 8./(25.*h2m);
//double b10 = sqrt(arg); //b^5 in units of sigma and epsilon
const double xStart = 0.7;

double delta_l[Lmax];
double Ea[2][sizeEa];
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
  return (-4/h2m*(x12 - x6) - l*(l+1)*x2 + Ex/h2m);//
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
  //parameters to determine change in sign of solution
  double k,x;
  double K;
  double j1[Lmax], j2[Lmax], n1[Lmax], n2[Lmax];
  double r1,r2;
  double R1[Lmax],R2[Lmax];
  double tmp0;
  char filename[17];
  int in=0;
  double sum=0;
  double pi = 4.*atan(1.);
  int fee=0;
  Emax = 3.5/epsilon +Estep;
  double Estart = 0.25/epsilon;
  int ee=0;
  for(E = Estart; E<Emax; E += Estep)
  {
    ee++;
    //printf("ee = %i\n", ee);
  }

  double sigmaTOT[ee];
  for(E = Estart; E<Emax; E += Estep)
  {
    k = sqrt(E/h2m);
    //scan on angular momentum
    for(l = 0; l<Lmax; l++)
    {
      tmp0 = integrate(E,l);
      double norm = CalcNorm();
      for(int j = 0; j < N; j++)
      {
        psi_norm[1][j] = psi[1][j]/norm;
        psi_norm[0][j] = psi[0][j];

        double RR1 = 4;
        //printf("%lg\n", RR1);

        if ((psi[0][j] > RR1)&&(psi[0][j]< RR1+3*h))
        {
          // while (fabs(psi_norm[1][j]) < 0.0001) {
          //   j++;
          // }
          do
          {
            r1 = psi[0][j];
            x = k*r1;
            Bessel(x,Lmax,0,j1);
            Bessel(x,Lmax,1,n1);
            R1[l] = psi[1][j];
            j++;
          } while((fabs(j1[l])<0.0005) || (fabs(n1[l])<0.0005));

          //printf("R1 = %g\n", R1[l]);
        }
        if ((psi[0][j] > r1)&&(psi[0][j]<r1+pi/(2*k)))
        {
          // while (fabs(psi_norm[1][j]) < 0.0001) {
          //   j++;
          // }

          do
          {
            r2 = psi[0][j];
            x = k*r1;
            Bessel(x,Lmax,0,j2);
            Bessel(x,Lmax,1,n2);
            R2[l] = psi[1][j];
            j++;
          } while((fabs(j2[l])<0.0005) || (fabs(n2[l])<0.0005));

        }
      }
      // printf("r1=%g    r2=%g\n", r1,r2);
      if(fee<11)
      {
        // sprintf(filename,"dati/plot%02i_%02i",l,fee);
        // writeCSVdouble(filename,(double *)psi_norm, 2, N);
      }

    }
    fee++;
    //calculate bessel function for given E
    x = k*r1;
    Bessel(x,Lmax,0,j1);
    Bessel(x,Lmax,1,n1);
    x = k*r2;
    Bessel(x,Lmax,0,j2);
    Bessel(x,Lmax,1,n2);
    //calculate delta_l

    for(l = 0; l<Lmax; l++)
    {
      // printf("l = %i\n", l);
       K = (R2[l]*r1)/(R1[l]*r2);
      //K = (R1[l])/(R2[l]);
      printf("fee = %i    K = %g\n",fee, K);
      delta_l[l] = atan2((K*j1[l] - j2[l]),(K*n1[l] - n2[l]));
      // printf("j1 = %g      j2 = %g      n1 = %g      n2 = %g\n",j1[l], j2[l], n1[l], n2[l] );
      // printf("delta_l = %g\n", delta_l[l]);
      sum = sum + (2*l+1)*sin(delta_l[l])*sin(delta_l[l]);
    }
    sigmaTOT[in] = 4.*pi/(k*k)*sum;
    Ea[0][in] = E*epsilon;
    Ea[1][in] = sigmaTOT[in];
    //printf("E = %g,    sigmaTOT_%i = %g\n",E, in, sigmaTOT[in]);
    sum = 0;
    in++;
  }

  //writeCSVdouble("CrossSec.csv",(double *)out, 2, ee);
  writeCSVdouble("energies.csv",(double *)Ea, 2, sizeEa);
  return 0;
}
