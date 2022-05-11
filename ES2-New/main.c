#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double pi=4.*atan(1.);
//# of basis functions
#define n 4
//Atomic Number
int Z = 1;
#define N 1
//gaussians exponents
double alpha1sH[4] = {13.00773, 1.962079, 0.444529, 0.1219492}; //values for Hydrogen
double alpha1sHe[4] = {14.899983, 2.726485, 0.757447, 0.251390};
double alpha1sBe[4] = {70.64859542, 12.92782254, 3.591490662, 1.191983464};
double alpha2sBe[4] = {3.072833610, 0.6652025433, 0.2162825386, 0.08306680972};
double alpha[n];

//matrixName
char s={'S'};char h={'H'};char c_n={'C'};char c_o={'O'};char f={'F'}; char de={'D'};
//mixing parameter
double mixing = 0.1;

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
void H(double Hm[n][n])
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
      Hm[p][q]= T_pq + A_pq;
    }
  }
  return;
}

//Direct and Exchange terms
void D(double Dm[n][n][n][n])
{
  double pi5m = sqrt(pi*pi*pi*pi*pi);
  double tmp1,tmp2;
  for(int p=0;p<n;p++)
  {
    for(int r=0;r<n;r++)
    {
      for(int q=0;q<n;q++)
      {
        tmp1 = alpha[p]+alpha[q];
        for(int s=0;s<n;s++)
        {
          tmp2 = alpha[r]+alpha[s];
          Dm[p][r][q][s] = 2*pi5m/(tmp1*tmp2*sqrt(tmp1+tmp2));
        }
      }
    }
  }
  return;
}
void FockMatrix(gsl_matrix *Fm, double Cm[N][n], double Dm[n][n][n][n])
{
  double F_pq = 0;
  for(int p=0;p<n;p++)
  {
    for(int q=0;q<n;q++)
    {
      F_pq = 0;
      for(int k=0;k<ceil(Z/2);k++)
      {
        for(int r=0;r<n;r++)
        {
          for(int s=0;s<n;s++)
          {
            F_pq += Cm[k][r]*Cm[k][s]*(2*Dm[p][r][q][s] - Dm[p][r][s][q]);
          }
        }
      }
    }
  }
}


//Definition of alpha based on the kind of Atom
void defAlpha()
{
  if (Z==1) for (int i = 0; i < n; i++) alpha[i] = alpha1sH[i];
  else if (Z==2) for (int i = 0; i < n; i++) alpha[i] = alpha1sHe[i];
  else if (Z==3)
  {
    for (int i = 0; i < n/2; i++) {
      alpha[i] = alpha1sBe[i];
      alpha[i+4] = alpha2sBe[i];
    }
  }
}
void printGslMatrix(char matrixName, gsl_matrix *m)
{
  int i,j;
  double evec[n][n];
  printf("%c matrix: \n", matrixName);
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      evec[i][j] = gsl_matrix_get(m,i,j);
      printf("%f  ", evec[i][j]);
    }
    printf("\n");
  }
}
void printDoubleMatrix(char matrixName, double M[n][n])
{
  int i,j;
  printf("%c matrix: \n", matrixName);
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      printf("%f  ", M[i][j]);
    }
    printf("\n");
  }
}
void printDMatrix(char matrixName, double M[n][n][n][n])
{
  int i,j;
  printf("%c matrix: \n", de);
  for(int p=0;p<n;p++)
  {
    printf("p: %d\n", p);
    for(int r=0;r<n;r++)
    {
      printf("r: %d\n", r);
      for(int q=0;q<n;q++)
      {
        for(int s=0;s<n;s++)
        {
          printf("%f ", Dm[p][r][q][s]);
        }
        printf("\n");
      }
    }
  }
  return;
}

//Real Generalized Symmetric-Definite Eigensystem
void RGSDE(gsl_matrix *Fm, gsl_matrix *Sm, gsl_matrix *Ev, gsl_vector *Ee)
{
  gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(4);
  gsl_eigen_gensymmv(Fm, Sm, Ee, Ev, w);
}

int main()
{
  //define alpha (gaussian exponents)
  defAlpha();
  //Matrices init
  double Hm[n][n];
  gsl_matrix *Sm = gsl_matrix_alloc(4,4);
  double Dm[n][n][n][n];
  double Cm[N][n];

  gsl_matrix *Ev = gsl_matrix_alloc(4,4);
  gsl_vector *Ee = gsl_vector_alloc(4);
  gsl_matrix_set_zero(Ev);
  gsl_vector_set_zero(Ee);

  S(Sm);
  printGslMatrix(s,Sm);
  H(Hm);
  printDoubleMatrix(h,Hm);
  D(Dm);

/*
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
  */
}
