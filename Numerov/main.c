
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"

#define N 200
#define Lmax 10
#define Nmax 10

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
  double x2 = x*x;
  return (-x2/4-l*(l+1)/x2+Ex);
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

  return (a[N-1]-a[N-3]);//
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
      for(int j = 0; j<4; j++)
      {
        double g = integrate(E);
        double g1 = integrate(E1);
        E2 = E;
        E = (E1*g-E*g1)/(g - g1);
        E1 = E2;

        if(abs(g) < 0)
        {
          printf("%i\n",j);
          break;
        }

        if(j>10)
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
