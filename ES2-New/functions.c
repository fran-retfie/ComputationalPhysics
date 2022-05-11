#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double pi=4.*atan(1.);
//# of basis functions
int n=4;
//Atomic Number
int Z = 1;
//gaussians exponents
double alpha1s[4] = {13.00773, 1.962079, 0.444529, 0.1219492}; //values for Hydrogen
//mixing parameter
double alpha = 0.1;

//overlap matrix (already integrated for GTO)
void S(gsl_matrix *Sm)
{
  int p,q;
  double tmp,tmp2;
  for(p=0;p<n;p++)
  {
    for(q=0;q<n;q++)
    {
      tmp2 = pi/(alpha[p] + alpha[q]);
      tmp = sqrt(tmp2*tmp2*tmp2);
      gsl_matrix_set(Sm,p,q,tmp);
    }
  }
  return;
}

//Hamiltonian matrix
void H(gsl_matrix *Hm)
{
  int p,q;
  double tmp,tmp5,T_pq,A_pq;
  double pi3 = pi*pi*pi;

  for(p=0;p<n;p++)
  {
    for(q=0;q<n;q++)
    {
      tmp = alpha[p] + alpha[q];
      tmp5 = tmp*tmp*tmp*tmp*tmp;
      T_pq = 3.*alpha[p]*alpha[q]*sqrt(pi3/tmp5);//Kin Energy
      A_pq = -2.*Z*pi/tmp; //Coulomb potential
      gsl_matrix_set(Hm,p,q,T_pq + A_pq);
    }
  }
  return;
}

//Direct and Exchange terms
void G(gsl_matrix *Gm)
{

}
