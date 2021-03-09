
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 500
#define Lmax 1
#define Nmax 2

const double tEnd = 5;
const double h = tEnd/(double) N;
const double Estep = 0.001;
double psi[4][N];
double Ea[2][Nmax*Lmax];
double a[N];
double b[N];

int l = 0;
double E = 2;

//effective potential (normalized!!)
double fr(double x, double Ex)
{
  if(x > 0) return (l*(l-1)/(x*x)-(x*x)+Ex);
  else return Ex;
}

double integrate(double Ex)
{
  double (*f)(double,double) = &fr;
  double gamma1;

  int sign = fmod(l,2) ? -1 : +1;

  //initial condition
  psi[0][0] = 0;
  psi[0][1] = h;
  psi[0][N-1] = h*(N-1);
  psi[0][N-2] = h*(N-2);
  a[0] = 0;
  a[1] = h;
  b[N-1] = sign*exp(-(N-1)*(N-1)*h*h/2)*psi[0][N-1];
  b[N-2] = sign*exp(-(N-2)*(N-2)*h*h/2)*psi[0][N-2];

  for(int j = 2; j < N; j++)
  {
    psi[0][j] = h*j;
    psi[0][N-j-1] = h*(N-j-1);
    a[j] = NumerovInt(psi[0][j], a[j-1], a[j-2], h, Ex, f);
    b[N-j-1] = NumerovIntBackward(psi[0][N-j-1], b[N-j], b[N-j+1], h, Ex, f);
    psi[1][j] = a[j]/psi[0][j];
    psi[2][N-j-1] = b[N-j-1]/psi[0][N-j-1];

    psi[3][j] = exp(-psi[0][j]*psi[0][j]/2);

    if(j>(N/2+1))
    {
      if((a[j-1]-b[j-1])*(a[j]-b[j])<0)
      {
        gamma1 = (a[j+1]-a[j-1]-b[j+1]+b[j-1])/psi[0][j];
        //break;
      }
      if((a[N-j-2]-b[N-j-2])*(a[N-j-1]-b[N-j-1])<0)
      {
        gamma1 = (a[N-j]-a[N-j-2]-b[N-j]+b[N-j-2])/psi[0][N-j-1];
        //break;
      }
    }
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];
  psi[2][N-2] = b[N-2]/psi[0][N-2];
  psi[2][N-1] = b[N-1]/psi[0][N-1];

  return gamma1;
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
    E = 2;
    //scan on energy
    for(int n = 0; n<Nmax; n++)
    {
      Eold = E;

      //simple scan
      for(; (gamma*gammaOld)>0|(E < Eold+Estep*3); E += Estep)
      {
        gammaOld = gamma;
        gamma = integrate(E);

        if(E>50)
        {
          printf("simple scan not converging, skipping!\n");
          break;
        }
      }

      //tangent method
      double E1 = E - Estep;
      double E2;
      for(int j = 0; fabs(E-E1)>Estep*0.1; j++)
      {
        double g = integrate(E);
        double g1 = integrate(E1);
        E2 = E;
        E = (E1*g-E*g1)/(g - g1);
        printf("%lg, %lg, %lg\n", E, g, g1);
        E1 = E2;

        if(j>4)
        {
          printf("tangent method not converging, skipping!\n");
          break;
        }
      }
      //if(psi[1][N-1] < 1)
      //{
        printf("%i, %i, %lg\n",l,n,E);
        Ea[1][n+l*Nmax] = E;
        Ea[0][n+l*Nmax] = (double) l;
      //}
      //else printf("skipped\n");
    }
  }


  writeCSVdouble("dati.csv",(double *)psi, 4, N);
  writeCSVdouble("energies.csv",(double *)Ea, 2, Nmax*Lmax);
  return 0;
}
