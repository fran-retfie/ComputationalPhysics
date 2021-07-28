#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double pi=4.*atan(1.);
int n=4;//# of basis functions
int Z=2;//# of electrons
//gaussians exponents
double alpha[4] = {14.899983, 2.726485, 0.757447, 0.251390};

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
      A_pq = -4.*pi/tmp; //Coulomb potential
      gsl_matrix_set(Hm,p,q,T_pq + A_pq);
    }
  }
  return;
}

//Direct matrix element
double D_pq(int p, int q, int r, int s)
{
  double pi5 = pi*pi*pi*pi*pi;
  double tmp = alpha[p] + alpha[q] + alpha[r] + alpha[s];
  double V_pqrs = 2.*sqrt(pi5)/((alpha[p] + alpha[s])*(alpha[r] + alpha[q])*sqrt(tmp));

  return V_pqrs;
}

//Exchange matrix element
double E_pq(int p, int q, int r, int s)
{
  double pi5 = pi*pi*pi*pi*pi;
  double tmp = alpha[p] + alpha[q] + alpha[r] + alpha[s];
  double V_pqrs = 2.*sqrt(pi5)/((alpha[p] + alpha[q])*(alpha[r] + alpha[s])*sqrt(tmp));

  return V_pqrs;
}

//Fock matrix
void F(gsl_matrix *Fm, gsl_matrix *Hm, gsl_matrix *Cm)
{
  int p,q,r,s,k;
  double H_pq, C_kr, C_ks;
  double F_pq = 0;
  double P_rs = 0;

  for(p=0;p<n;p++)
  {
    // printf("p=%d\n", p);
    for(q=0;q<=p;q++)
    {
      // printf("q=%d\n", q);
      for(r=0;r<p;r++)
      {
        // printf("r=%d\n", r);
        if(r<p-1)
        {
          for(s=0;s<r;s++)
          {
            P_rs = 0;
            // printf("s=%d\n", s);
            for(k=0;k<=Z/2;k++)
            {
              C_kr = gsl_matrix_get(Cm,k,r);
              C_ks = gsl_matrix_get(Cm,k,s);
              P_rs = P_rs + C_kr*C_ks;//NOTE: sum over k up to Z/2
            }
            // printf("%f\n", P_rs);
            F_pq = F_pq + P_rs*(D_pq(p,q,r,s) + E_pq(p,q,r,s));
          }
        }
        else if(r==p-1)
        {
          for(s=0;s<q;s++)
          {
            P_rs = 0;
            // printf("s=%d\n", s);
            for(k=0;k<=Z/2;k++)
            {
              C_kr = gsl_matrix_get(Cm,k,r);
              C_ks = gsl_matrix_get(Cm,k,s);
              P_rs = P_rs + C_kr*C_ks;//sum over k up to N/2 (NOTA: quindi le ultime due righe dei coefficienti non servono a niente?)
            }
            // printf("%f\n", P_rs);
            F_pq = F_pq + P_rs*(D_pq(p,q,r,s) + E_pq(p,q,r,s));
          }
        }
      }
      H_pq = gsl_matrix_get(Hm,p,q);
      F_pq = F_pq + H_pq;
      // printf("p=%d, q=%d, F_pq=%f\n",p,q, F_pq);
      gsl_matrix_set(Fm,p,q,F_pq);
      gsl_matrix_set(Fm,q,p,F_pq);
    }
  }
  return;
}

//Initial coefficient guess (such that sum_{k=0 to Z/2})
// void ICG(gsl_matrix *Cm)
// {
//   // gsl_matrix_set(Cm,0,0,-1);gsl_matrix_set(Cm,0,1,-1);gsl_matrix_set(Cm,0,2,1);gsl_matrix_set(Cm,0,3,-1);
//   // gsl_matrix_set(Cm,1,0,-1);gsl_matrix_set(Cm,1,1,1);gsl_matrix_set(Cm,1,2,1);gsl_matrix_set(Cm,1,3,1);
//   // gsl_matrix_set(Cm,2,0,1);gsl_matrix_set(Cm,2,1,1);gsl_matrix_set(Cm,2,2,-1);gsl_matrix_set(Cm,2,3,1);
//   // gsl_matrix_set(Cm,3,0,1);gsl_matrix_set(Cm,3,1,-1);gsl_matrix_set(Cm,3,2,-1);gsl_matrix_set(Cm,3,3,-1);
//   gsl_matrix_set(Cm,0,0,0);gsl_matrix_set(Cm,0,1,0);gsl_matrix_set(Cm,0,2,0);gsl_matrix_set(Cm,0,3,0);
//   gsl_matrix_set(Cm,1,0,0);gsl_matrix_set(Cm,1,1,0);gsl_matrix_set(Cm,1,2,0);gsl_matrix_set(Cm,1,3,0);
//   gsl_matrix_set(Cm,2,0,0);gsl_matrix_set(Cm,2,1,0);gsl_matrix_set(Cm,2,2,0);gsl_matrix_set(Cm,2,3,0);
//   gsl_matrix_set(Cm,3,0,0);gsl_matrix_set(Cm,3,1,0);gsl_matrix_set(Cm,3,2,0);gsl_matrix_set(Cm,3,3,0);
// }

double minE(gsl_vector *Ee)
{
  int i;
  double tmp,min;
  min = gsl_vector_get(Ee,0);
  for(i=1;i<n;i++)
  {
    tmp = gsl_vector_get(Ee,i);
    if(tmp<min)
    {
      min = tmp;
    }
  }
  return min;
}

//Real Generalized Symmetric-Definite Eigensystem
void RGSDE(gsl_matrix *Fm, gsl_matrix *Sm, gsl_matrix *Cm, gsl_vector *Ee)
{
  gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);
  gsl_eigen_gensymmv(Fm, Sm, Ee, Cm, w);
  gsl_eigen_gensymmv_free(w);
  return;
}

int main()
{
  gsl_matrix *Hm = gsl_matrix_alloc(n,n);
  gsl_matrix *Sm = gsl_matrix_alloc(n,n);
  gsl_matrix *Fm = gsl_matrix_alloc(n,n);
  gsl_matrix *C_old = gsl_matrix_alloc(n,n);
  gsl_matrix *C_new = gsl_matrix_alloc(n,n);
  gsl_vector *Ee = gsl_vector_alloc(n);
  gsl_matrix_set_zero(Hm);
  gsl_matrix_set_zero(Sm);
  gsl_matrix_set_zero(Fm);
  gsl_matrix_set_zero(C_old);
  gsl_matrix_set_zero(C_new);
  gsl_vector_set_zero(Ee);

  S(Sm); //build overlap matrix
  H(Hm); //build Hamiltonian

  //self-consistent procedure
  double alfa = 1e-3;//coefficient in the variation of the coefficients
  double check = 1;
  double convergence = 1e-6;
  double tmp;
  int i,j;
  //first iteration (to set the min of energy)
  F(Fm,Hm,C_new); //build Fock matrix
  RGSDE(Hm,Sm,C_new,Ee); //Solve the generalized eigrnvalue problem
  check = abs(minE(Ee));//POSSIBLE PROBLEM IF ALREADY SMALLER THAN convergence
  printf("check = %f\n", check);
  while(check>convergence)
  {
    //compute new coefficients as a small variation of
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        tmp = alfa*gsl_matrix_get(C_new,i,j) + (1-alfa)*gsl_matrix_get(C_old,i,j);
        gsl_matrix_set(C_new,i,j,tmp);
      }
    }
    C_old = C_new;
    F(Fm,Hm,C_new); //build Fock matrix
    RGSDE(Hm,Sm,C_new,Ee); //Solve the generalized eigrnvalue problem
    tmp = minE(Ee);
    check = abs(check-tmp);
    printf("check = %f\n", check);
  }




  double evec[n][n], eval[n];

  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      evec[i][j] = gsl_matrix_get(C_new,i,j);
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
