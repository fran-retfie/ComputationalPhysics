#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

const double pi=4.*atan(1.);
//# of basis functions
#define n 8
//Atomic Number
int Z = 4;
//Number of electrons
#define N 4
//Number of orbitals
int orbitalNum = 2;
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
//Direct term
void D(double Dm[n][n][n][n])
{
  double pi5m = sqrt(pi*pi*pi*pi*pi);
  double tmp1,tmp2;
  for(int p=0;p<n;p++)
  {
    for(int q=0;q<n;q++)
    {
      for(int r=0;r<n;r++)
      {
        tmp1 = alpha[p]+alpha[q];
        for(int s=0;s<n;s++)
        {
          tmp2 = alpha[r]+alpha[s];
          Dm[p][q][r][s] = 2*pi5m/(tmp1*tmp2*sqrt(tmp1+tmp2));
        }
      }
    }
  }
  return;
}
//Exchange term
void Ex(double Exm[n][n][n][n])
{
  double pi5m = sqrt(pi*pi*pi*pi*pi);
  double tmp1,tmp2;
  for(int p=0;p<n;p++)
  {
    for(int q=0;q<n;q++)
    {
      for(int r=0;r<n;r++)
      {
        tmp1 = alpha[p]+alpha[r];
        for(int s=0;s<n;s++)
        {
          tmp2 = alpha[q]+alpha[s];
          Exm[p][q][r][s] = 2*pi5m/(tmp1*tmp2*sqrt(tmp1+tmp2));
        }
      }
    }
  }
  return;
}
void FockMatrix(gsl_matrix *Fm, gsl_matrix *Cm, double Exm[n][n][n][n], double Dm[n][n][n][n], double Hm[n][n])
{
  double F_pq = 0;
  for(int p=0;p<n;p++)
  {
    for(int q=0;q<n;q++)
    {
      F_pq = 0;
      for(int k=0;k<orbitalNum;k++)
      {
        for(int r=0;r<n;r++)
        {
          for(int s=0;s<n;s++)
          {
            F_pq += gsl_matrix_get(Cm, k,r)*gsl_matrix_get(Cm, k,s)*(2*Dm[p][q][r][s] - Exm[p][q][r][s]);
          }
        }
      }
      gsl_matrix_set(Fm,p,q,Hm[p][q]+F_pq);
    }
  }
}
//Definition of alpha based on the kind of Atom
void defAlpha()
{
  if (Z==1) for (int i = 0; i < n; i++) alpha[i] = alpha1sH[i];
  else if (Z==2) for (int i = 0; i < n; i++) alpha[i] = alpha1sHe[i];
  else if (Z==4)
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
          printf("%f ", M[p][r][q][s]);
        }
        printf("\n");
      }
    }
  }
  return;
}
void printCMatrix(double M[N][n])
{
  printf("Cm matrix: \n");
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<n;j++)
    {
      printf("%f  ", M[i][j]);
    }
    printf("\n");
  }
}
//Hartree-Fock Energy
void E_HF(gsl_matrix Cm, double eval[orbitalNum], double Dm[n][n][n][n], double Exm[n][n][n][n], double Hm[n][n], gsl_matrix *E)
{
  double result=0;
  double tmp = 0;
  for(int i=0;i<orbitalNum;i++)
  {
    for(int k=0;k<orbitalNum;k++)
    {
      for(int p=0;p<n;p++)
      {
        for(int q=0;q<n;q++)
        {
          for(int r=0;r<n;r++)
          {
            for(int s=0;s<n;s++)
            {
              tmp += gsl_matrix_get(Cm, i,p)*gsl_matrix_get(Cm, k,r)*gsl_matrix_get(Cm, i,q)*gsl_matrix_get(Cm, k,s)*(2*Dm[p][q][r][s]-Exm[p][q][r][s]); //No 2 for H
            }
          }
        }
      }
    }
  }
  double tmp2 = 0;
  for(int i=0;i<orbitalNum;i++)
  {
    for(int p=0;p<n;p++)
    {
      for(int q=0;q<n;q++)
      {
        tmp2 += Cm[i][p]*Cm[i][q]*Hm[p][q];
      }
    }
  }


  for(int r=0;r<orbitalNum;r++)
  {
    result += eval[r];
  }
  E[0] = 2*result - tmp;//result-tmp/2; for H
  E[1] = 2*tmp2 + tmp;//tmp2 + tmp/2; for H (not working)
}
/* generate a random floating point number from min to max */
double randfrom(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
void initCoef(double Cm[N][n])
{
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<n;j++)
    {
      Cm[i][j] = randfrom(0.0,1.0);
    }
  }
}
void updateCoef(gsl_matrix *Cm, gsl_matrix *CmNew)
{
  double tmp=0;
  for(int k=0;k<orbitalNum;k++)//column
  {
    for(int i=0;i<n;i++)//row
    {
      tmp=mixing*gsl_matrix_get(CmNew,i,k) + (1-mixing)*gsl_matrix_get(Cm,i,k);
      gsl_matrix_set(Cm,i,k,tmp);
    }
  }
}
//tks volpx
void Roothan(gsl_matrix *C, gsl_matrix *E, const gsl_matrix *F,const gsl_matrix *V,gsl_eigen_symmv_workspace *ws, gsl_matrix *A, gsl_matrix *Fp)
{
  gsl_matrix_set_zero(A);
  gsl_matrix_set_zero(Fp);

  //Compute Fp=F'=V^T F V
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, F, 0.0, A);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, V, 0.0, Fp);

  //Sort eigenvalue problem for Fp
  gsl_eigen_symmv(Fp, E, A, w);
  gsl_eigen_symmv_sort(E,A,GSL_EIGEN_SORT_VAL_ASC);

  //C=V C'
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, A, 0.0, C);

  return;
}



int main()
{
  //seed
  srand (time ( NULL));
  //define alpha (gaussian exponents)
  defAlpha();
  for(int i=0;i<n;i++)
  {
    printf("%lf\n", alpha[i]);
  }
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  //Matrices init
  double Hm[n][n];
  gsl_matrix *Sm = gsl_matrix_alloc(n,n);
  gsl_matrix *SVec = gsl_matrix_alloc(n,n);
  gsl_matrix *SVal = gsl_vector_alloc(n);
  double Dm[n][n][n][n];
  double Exm[n][n][n][n];
  gsl_matrix *Cm = gsl_matrix_alloc(n,n);
  gsl_matrix *CmNew = gsl_matrix_alloc(n,n);
  double Eval[orbitalNum];

  gsl_matrix *X = gsl_matrix_alloc(n,n);
  gsl_matrix *E = gsl_matrix_alloc(n,n);
  gsl_matrix *Fm = gsl_matrix_alloc(n,n);
  gsl_matrix *evec = gsl_matrix_alloc(n,n);
  gsl_vector *eval = gsl_vector_alloc(n);
  gsl_matrix_set_zero(Sm);
  gsl_matrix_set_zero(SVec);
  gsl_vector_set_zero(SVal);
  gsl_matrix_set_zero(Cm);
  gsl_matrix_set_zero(CmNew);
  gsl_matrix_set_zero(X);
  gsl_matrix_set_zero(E);
  gsl_matrix_set_zero(evec);
  gsl_matrix_set_zero(Fm);
  gsl_vector_set_zero(eval);

  double thr = 0.00001;
  int counter = 0;

  H(Hm);
  S(Sm);
  gsl_eigen_symmv(Sm,SVal,SVec,w);
  gsl_eigen_symmv_sort(SVal,SVec,GSL_EIGEN_SORT_VAL_ASC);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      gsl_matrix_set(X,i,j,gsl_matrix_get(SVec,i,j) / gsl_matrix_get(SVal,i,j));
    }
  }
  D(Dm);
  Ex(Exm);

  printDoubleMatrix(h,Hm);
  printGslMatrix(s,Sm);
  //printDMatrix(de,Dm);

  double tmp,tmpOld,Eold,Enew;
  int index[orbitalNum];
  Eold=1;
  Enew=0;

  while(fabs(Enew-Eold)>thr)
  //for(int ii=0;ii<200;ii++)
  {
    counter++;
    FockMatrix(Fm,Cm,Dm,Hm);
    Roothan(Cm,E,Fm,X,w,A,Fp);
    }
    //Get new coeff from evec

    updateCoef(Cm,CmNew);
    Eold = gal_matrix_get(E,0,0);
    E_HF(Cm,Eval,Dm,Hm,E);
    printf("E1: %lf     E2: %lf\n", E[0],E[1]);


  }
  gsl_eigen_symmv_free(w);
  printf("E: %lf    iterations: %d\n", E[0], counter);
}
