
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>

#include "CSVio.h"
#include "DiffInt.h"
#include "doubleDef.h"

//------------------------------------------------------------------------------

//grid size
#define N 1000
//max number of shells
#define Smax 6
//max of quantum number n
#define Nmax 2
//max of quantum number l
#define Lmax 4
//max number of iteration
#define Qmax 150

const int Nelist[Smax] = {2,8,18,20,26,40};
const int nnList[Smax] = {0,1*Nmax,2*Nmax,1,1+1*Nmax,3*Nmax};

const double_t tEnd = 27;
const double_t tStart = tEnd/(double_t) N;
const double_t h = (tEnd-tStart)/(double_t) N;
const double_t Estep = 0.003;
const double_t h2_12 = h*h/12;

double_t Vext[N];
double_t psi[3][N];
double_t psimem[Nmax*Lmax][N];
double_t rho[4][N];
double_t rhoOld[N];
double_t Ea[3][Nmax*Lmax];
//bool Eavail[Nmax*Lmax];

char title[30];
char filename[30];

int q = 0;
int n = 1;
int l = 0;
double_t E = 0;

int ne;
int nn;

//state mixing coefficient
const double_t mix = 0.15;

//correlation terms constants
const double_t A = 0.031091;
const double_t alpha1 = 0.21370;
const double_t beta1 = 7.5957;
const double_t beta2 = 3.5876;
const double_t beta3 = 1.6382;
const double_t beta4 = 0.49294;

//static const double_t coeff3 = 3/8*h;
static const double_t coeff2 = h/3;

//Na atom
const double_t rs = 3.93;
//K atom
//const double_t rs = 4.86;

int Ne;
//const double_t rs3 = (rs*rs*rs);
//const double_t o_rs3 = 1/rs3;
double_t Rc;
double_t Rc2;
double_t Rc3;

double_t C1;
double_t C2;
const double_t C3 = 2*A/3;

double_t gammaNew;
double_t gammaOld;
double_t Eold;

//------------------------------------------------------------------------------

double_t get_r(int k)
{
  return tStart + h*k;
}

double_t get_r_max(int j, int i)
{
  return ( (i>j) ? (tStart + h*i) : (tStart + h*j) );
}

//calculate the Hartree term
double_t CalcHartree(int nr)
{
  double_t I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += (get_r(j-1)*get_r(j-1)*rho[1][j-1]/get_r_max(j-1,nr) + 4*get_r(j)*get_r(j)*rho[1][j]/get_r_max(j,nr) + get_r(j+1)*get_r(j+1)*rho[1][j+1]/get_r_max(j+1,nr) );
  }

  return 4*M_PI * h/3 * I;
}

//effective potential (normalized!!)
double_t fr(int k)
{
  double_t r = get_r(k);
  double_t r2 = r*r;

  double_t rse = cbrtq(3/(4*M_PI*rho[1][k]));
  double_t rse2 = rse*rse;
  double_t rse05 = sqrtq(rse);
  double_t rse15 = rse05*rse;

  double_t rsenum = 0.5*beta1*rse05 + beta2*rse + 1.5*beta3*rse15 + 2*beta4*rse2;
  double_t rsefac = beta1*rse05 + beta2*rse + beta3*rse15 + beta4*rse2;

  double_t Vext =  ( (r < Rc) ? (C1*(r2 - 3*Rc2)) : (-Ne/r) );
  double_t Vhart = (q != 0) ? CalcHartree(k) : 0;
  double_t Vexc =  (q != 0) ? C2*cbrtq(rho[1][k]) : 0;
  double_t Vcorr = (q != 0) ? -C3*( (2*alpha1*rs+3)*logq(1+1/(2*A*rsefac)) + (alpha1*rs+1)*rsenum/((2*A*rsefac+1)*rsefac) ): 0;

  return ( -l*(l+1)/r2 - 2*( Vhart + Vexc + Vcorr + Vext));
}

//integrate Schrodinger equation with Numerov
double_t integrate(double_t Ex)
{
  psi[1][0] = h;
  psi[1][1] = 2*h;

  for(int j = 2; j < N; j++) psi[1][j] = ((2-10*(Vext[j-1]+Ex)*h2_12)*psi[1][j-1] - (1+(Vext[j-2]+Ex)*h2_12)*psi[1][j-2])/(1+(Vext[j]+Ex)*h2_12);

  return (psi[1][N-1] - psi[1][N-2]);
}

//calculate the norm of wavefunction
double_t CalcNorm()
{
  double_t I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += (get_r(j-1)*get_r(j-1)*psi[1][j-1]*psi[1][j-1] + 4*get_r(j)*get_r(j)*psi[1][j]*psi[1][j] + get_r(j+1)*get_r(j+1)*psi[1][j+1]*psi[1][j+1]);
  }
  return sqrtq(h/3*I);
}

//------------------------------------------------------------------------------

int main()
{
  for(int j = 0; j < N; j++) psi[0][j] = get_r(j);
  for (int k = 0; k < N; k++) rho[1][k] = 0;

  //loop over the first Nmax sheels
  for(int s = 1; s <= Smax; s++)
  {
    //calculate the number of electrons and Rc
    Ne = Nelist[s-1];
    Rc = cbrtq(Ne)*rs;
    Rc2 = Rc*Rc;
    Rc3 = Rc*Rc*Rc;
    C1 = 0.5*Ne/Rc3;
    C2 = cbrtq(3/M_PI);

    //khon-sham loop
    for (q = 0; q < Qmax; q++)
    {

      //scan on angular momentum
      for(l = 0; l < Lmax; l++)
      {
        //start research form E = -10
        E = -10;

        //scan on energy over (Ne/2)+1) states
        for(n = 1; n <= Nmax ; n++)
        {
          //calculate potential
          for (int i = 0; i < N; i++) Vext[i] = fr(i);

          //simple scan
          Eold = E;
          for(; (gammaNew*gammaOld)>0|(E < Eold+Estep*4); E += Estep)
          {
            gammaOld = gammaNew;
            gammaNew = integrate(E);
            if(E>0) break;
          }

          double_t E1 = E - 2*Estep;
          double_t E2;
          double_t g;
          double_t g1;

          //tangent method
          for(int j = 0; j<50; j++)
          {
            g = integrate(E);
            g1 = integrate(E1);
            //if(j>200) printf("%g\n",g);
            if(fabsq(g-g1)<1e-20) break;
            E2 = E;
            E = (E1*g-E*g1)/(g - g1);
            E1 = E2;
            //if(fabsl(g-g1)<1e-200) break;
          }

          /*
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
          */

          //printf("Ne = %i,\tl = %i,\tn = %i,\tE = %40.39lg\n",Ne, l, n, E);
          Ea[0][(n-1) + l*Nmax] = (double_t) n;
          Ea[1][(n-1) + l*Nmax] = (double_t) l;
          Ea[2][(n-1) + l*Nmax] = 2*E;
          //Eavail[(n-1) + l*Nmax] = false;

          //save the potential
          for (int i = 0; i < N; i++) psi[2][i] = fr(i);

          //find radial wavefunction
          for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/psi[0][j];

          double_t norm = CalcNorm() * 2*M_SQRTPI; //normalize to 1/sqrt(4*pi)
          for(int j = 0; j < N; j++)
          {
            psi[1][j] = psi[1][j]/norm;
            psimem[(n-1) + l*Nmax][j] = psi[1][j];
          }

          sprintf(filename, "dati/plot%02i%02i%02i", s, n, l);
          writeCSVdouble_t(filename, (double_t *) psi, 3, N, title);
        }
      }

      sprintf(filename, "dati/energies%02i", s);
      writeCSVdouble_t(filename, (double_t *) Ea, 3, Nmax*Lmax, title);

//---------------------------density stuff---------------------------------------

      //now calculate the density here
      //search for lowest energy states
      ne = 0;
      for (int k = 0; k < N; k++) rhoOld[k] = rho[1][k];
      for (int k = 0; k < N; k++) rho[1][k] = rho[1][k] * (1 - mix);
      nn = 0;

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

        for(int k = 0; k < N; k++)
        {
          double_t val = psimem[(nmin-1) + lmin*Nmax][k];
          rho[0][k] = psi[0][k];
          rho[1][k] += 2*(2*lmin + 1)*val*val * mix;
        }
      }

      /*
      for(int k = 0; k < N; k++)
      {
        rho[2][k] = CalcHartree(k); //hartree term
        rho[3][k] = cbrtq(3/(4*M_PI*rho[1][k])); //r_s
      }
      */
    }

    sprintf(title, "N_e = %i", ne);
    sprintf(filename, "dati/density%02i", s);
    writeCSVdouble_t(filename, (double_t *) rho, 4, N, title);
    printf("\n");
  }

  return 0;
}
