#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double pi=4.*atan(1.);
int n=4;//# of basis functions
//gaussians exponents
double alpha[4] = {13.00773, 1.962079, 0.444529, 0.1219492};

//gaussian basis
double Chi(double r, double alpha)
{
  double r2 = r*r;

  return (exp(-alpha*r2));
}

//overlap matrix (already integrated for GTO)
void S(gsl_matrix *Sm)
{
  int p,q;
  double tmp;
  for(p=0;p<n;p++)
  {
    for(q=0;q<n;q++)
    {
      tmp = sqrt((pi/(alpha[p] + alpha[q]))*(pi/(alpha[p] + alpha[q]))*(pi/(alpha[p] + alpha[q])));
      // printf("%f\n", tmp);
      gsl_matrix_set(Sm,p,q,tmp);
      // printf("gslMatrix: %f\n", gsl_matrix_get(Sm,p,q));
      // printf("-------------------------------------------\n");
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
      T_pq = 3.*alpha[p]*alpha[q]*sqrt(pi3/tmp5);
      A_pq = -2.*pi/tmp; //Coulomb potential
      gsl_matrix_set(Hm,p,q,T_pq + A_pq);
    }
  }
  return;
}

//Real Generalized Symmetric-Definite Eigensystem
void RGSDE(gsl_matrix *Hm, gsl_matrix *Sm, gsl_matrix *Ev, gsl_vector *Ee)
{
  gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(4);
  gsl_eigen_gensymmv(Hm, Sm, Ee, Ev, w);
}

int main()
{
  gsl_matrix *Hm = gsl_matrix_alloc(4,4);
  gsl_matrix *Sm = gsl_matrix_alloc(4,4);
  gsl_matrix *Ev = gsl_matrix_alloc(4,4);
  gsl_vector *Ee = gsl_vector_alloc(4);
  gsl_matrix_set_zero(Hm);
  gsl_matrix_set_zero(Sm);
  gsl_matrix_set_zero(Ev);
  gsl_vector_set_zero(Ee);

  S(Sm); //build overlap matrix
  H(Hm); //build Hamiltonian
  RGSDE(Hm,Sm,Ev,Ee); //Solve the generalized eigrnvalue problem

  double evec[4][4], eval[4];
  int i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      evec[i][j] = gsl_matrix_get(Ev,i,j);
      printf("%f  ", evec[i][j]);
    }
    printf("\n");
  }
  printf("----------------------------------\n");
  for(i=0;i<n;i++)
  {
    eval[i] = gsl_vector_get(Ee,i);
    printf("%f  ", eval[i]);
    printf("\n");
  }
}
