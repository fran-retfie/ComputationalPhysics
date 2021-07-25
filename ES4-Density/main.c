
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <stdlib.h>
#include "CSVio.h"
#include "DiffInt.h"
#include <stdbool.h>
#include "doubleDef.h"

//grid size
#define N 1000

//max number of shells
#define Smax 4

//max of quantum number n
#define Nmax 4

//max of quantum number l
#define Lmax 4

const int Nelist[Smax] = {2,8,18,20};
const int nnList[Smax] = {0,1*Nmax,2*Nmax,1};

const double_t tEnd = 25;
const double_t h = tEnd/(double_t) N;
const double_t Estep = 0.001;
double_t psi[3][N];
double_t a[N];
double_t psimem[Nmax*Lmax][N];
double_t rho[2][N];
double_t Ea[3][Nmax*Lmax];
bool Eavail[Nmax*Lmax];

char title[30] = "prova";

int n = 1;
int l = 0;
double_t E = 0;

//Na atom
const double_t rs = 3.93;
//K atom
//const double_t rs = 4.86;

const double_t rs3 = (rs*rs*rs);
const double_t o_rs3 = 1/rs3;
double_t Rc;
double_t Rc2;
double_t Rc3;

//effective potential (normalized!!)
double_t fr(double_t r, double_t Ex)
{
  double_t r2 = r*r;
  double_t Vext = o_rs3 * ( (r < Rc) ? (r2 - 3*Rc2) : (-2*Rc3/r) );
  return (-Vext -l*(l+1)/r2 + Ex);
}

//integrate Schrodinger equation with Numerov
double_t integrate(double_t Ex)
{
  double_t (*f)(double_t,double_t) = &fr;

  psi[0][0] = 0;
  a[0] = 0;
  psi[0][1] = h;
  a[1] = h;
  psi[0][2] = 2*h;
  a[2] = 2*h;

  for(int j = 3; j < N; j++)
  {
    psi[0][j] = h*j;
    a[j] = NumerovInt_t(psi[0][j-1], a[j-1], a[j-2], h, Ex, f);
    psi[1][j] = a[j];
  }

  psi[1][0] = psi[1][2];
  psi[1][1] = psi[1][2];

  return (a[N-1]-a[N-2]);
}

//calculate the norm of wavefunction
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
  //loop over the first Nmax sheels
  for(int s = 1; s <= Smax; s++)
  {
    //calculate the number of electrons and Rc

    int Ne = Nelist[s-1];
    Rc = cbrtq(Ne)*rs;
    Rc2 = Rc*Rc;
    Rc3 = Rc*Rc*Rc;

    double_t gamma = 1;
    double_t gammaOld = 1;
    double_t Eold;

    //scan on angular momentum
    for(l = 0; l < Lmax; l++)
    {
      //start research form E = -10
      E = -10;

      //scan on energy over (Ne/2)+1) states
      for(n = 1; n <= Nmax ; n++)
      {
        Eold = E;
        //simple scan
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

        //bisection method
        for(int j = 0; j<200; j++)
        {
          E2 = (E+E1)*0.5;
          g1 = integrate(E2);
          g = integrate(E);
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
          //if(fabs(g-g1)<1e-100) break;
          //if(j>200) printf("%g\n",E);
        }

        printf("Ne = %i,\tl = %i,\tn = %i,\tE = %40.39lg\n",Ne, l, n, E);
        Ea[0][(n-1) + l*Nmax] = (double_t) n;
        Ea[1][(n-1) + l*Nmax] = (double_t) l;
        Ea[2][(n-1) + l*Nmax] = E;
        Eavail[(n-1) + l*Nmax] = false;

        double_t norm = CalcNorm();
        for(int j = 0; j < N; j++)
        {
          psi[1][j] = psi[1][j]/norm;
          psimem[(n-1) + l*Nmax][j] = psi[1][j];
        }

        char filename[30];
        sprintf(filename, "dati/plot%02i%02i%02i", s, n, l);
        writeCSVdouble_t(filename, (double_t *) psi, 3, N, title);
      }
    }

    char filename[30];
    sprintf(filename, "dati/energies%02i", s);
    writeCSVdouble_t(filename, (double_t *) Ea, 3, Nmax*Lmax, title);

    //now calculate the density here
    //search for lowest energy states
    int ne = 0;
    for (int k = 0; k < N; k++) rho[1][k] = 0;

    int nn = 0;

    while(ne < Ne)
    {
      /*
      //find minimum energy
      int  jmin = Nmax*Lmax-1;
      for (int j = 0; j < Nmax*Lmax; j++)
        if((Ea[2][jmin] > Ea[2][j]) && !Eavail[j]) jmin = j;
      Eavail[jmin] = true;
      */

      int jmin = nnList[nn];
      nn++;

      int nmin = (int) Ea[0][jmin];
      int lmin = (int) Ea[1][jmin];
      ne += 2*(2*lmin + 1);
      printf("ne = %i,\tn = %i,\tl = %i,\tE = %lg\n", ne, nmin, lmin, Ea[2][jmin]);

      for (int k = 0; k < N; k++)
      {
        double_t val = psimem[(nmin-1) + lmin*Nmax][k];
        rho[0][k] = psi[0][k];
        rho[1][k] += 2*(2*lmin + 1)*val*val;
      }
    }

    sprintf(title, "N_e = %i", ne);
    sprintf(filename, "dati/density%02i", s);
    writeCSVdouble_t(filename, (double_t *) rho, 2, N, title);
    printf("\n");
  }
  return 0;
}
