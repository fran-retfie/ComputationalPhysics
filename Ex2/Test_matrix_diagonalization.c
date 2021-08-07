#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

double pi=4.*atan(1.);
int n=4;//# of basis functions
int Z=2;//# of electrons

double alpha[4] = {14.899983, 2.726485, 0.757447, 0.251390}; //gaussians exponents (prof)
// double alpha[4] = {38.216677, 5.749982, 1.236745, 0.297104}; //gaussians exponents (web) BEST
// double alpha[4] = {38.421634, 5.77803, 1.241774, 0.297964}; //gaussians exponents (web)
char s={'S'};char h={'H'};char c_n={'C'};char c_o={'O'};char f={'F'}; //used in printMatrix (matrix name)

void printMatrix(int matrixName, gsl_matrix *m)
{
  int i,j;
  double evec[n][n];
  printf("--------------------------- %d ---------------------------\n", matrixName);
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      evec[i][j] = gsl_matrix_get(m,i,j);
      printf("%f  ", evec[i][j]);
    }
    printf("\n");
  }
  printf("---------------------------  ---------------------------\n");
}
void printVector(int matrixName, gsl_vector *m)
{
  int i;
  double evec[n];
  printf("--------------------------- %d ---------------------------\n", matrixName);
  for(i=0;i<n;i++)
  {
    evec[i] = gsl_vector_get(m,i);
    printf("%f  ", evec[i]);
    printf("\n");
  }
  printf("---------------------------  ---------------------------\n");
}

void matrixProduct(gsl_matrix *A, gsl_matrix *B, gsl_matrix *Out)
{
  int i,j,k;
  double a_ij,b_ji,tmp;
  for(i=0;i<n;i++)//Out row
  {
    for(j=0;j<n;j++)//Out column
    {
      tmp = 0;
      for(k=0;k<n;k++)
      {
        a_ij = gsl_matrix_get(A,i,k);
        b_ji = gsl_matrix_get(B,k,j);
        tmp += a_ij*b_ji;
      }
      gsl_matrix_set(Out,i,j,tmp);
    }
  }
}

//overlap matrix (already integrated for GTO)
void S(gsl_matrix *Sm, gsl_matrix *evec)
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
  printMatrix(1,Sm);
  return;
}

void matrixDiag(gsl_matrix *A, gsl_matrix *Aevec, gsl_matrix *diagMat)
{
  int p,q;
  double tmp;
  //diagonalization
  gsl_vector *eval = gsl_vector_alloc(n);
  gsl_matrix *Acopy = gsl_matrix_alloc(n,n);
  gsl_matrix_memcpy(Acopy,A);
  gsl_vector_set_zero(eval);
  gsl_matrix_set_zero(Aevec);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(Acopy,eval,Aevec,w);

  // sort eigenvalues and eigenvectors for increasing eigenvalues
  double perm[n];
  gsl_vector *index = gsl_vector_alloc(n);
  gsl_vector *order = gsl_vector_alloc(n);
  for(p=0;p<n;p++)
  {
    gsl_vector_set(index,p,p);
  }
  gsl_sort_vector2(eval,index);
  double pq;int i;
  gsl_matrix *tmp3 = gsl_matrix_alloc(n,n);
  for(p=0;p<n;p++)
  {
    i = (int)gsl_vector_get(index,p);
    for(q=0;q<n;q++)
    {
      gsl_vector_set(order,q,gsl_matrix_get(Aevec,q,i));
    }
    gsl_matrix_set_col(tmp3,p,order);
  }

  // build the diagonal matrix
  // gsl_matrix_set_zero(A);
  tmp=0;
  double tmp2 = 0;
  for(p=0;p<n;p++)
  {
    tmp += gsl_vector_get(eval,p)*gsl_vector_get(eval,p);
    gsl_matrix_set(diagMat,p,p,gsl_vector_get(eval,p));
  }
  tmp = 1./sqrt(tmp);
  //normalization
  // gsl_matrix_scale(A,tmp);
  // gsl_matrix_scale(evec,tmp);
}

// diagonalizing matrix : X -> X^AX = 1
void diagonalizingMatrix(gsl_matrix *diagMat, gsl_matrix *evec, gsl_matrix *X)
{
  gsl_matrix *s_half = gsl_matrix_alloc(n,n);
  gsl_matrix *tmp = gsl_matrix_alloc(n,n);
  // Canonical Orthogonalization: X = Us^-1/2 (U = evec)
  for(int i=0;i<n;i++)
  {
    gsl_matrix_set(s_half,i,i,1./sqrt(gsl_matrix_get(diagMat,i,i)));
  }
  matrixProduct(evec,s_half,X);
  //normalize columns
  // double X_ji,sum;
  // for(int i=0;i<n;i++)
  // {
  //   sum = 0;
  //   for(int j=0;j<n;j++)
  //   {
  //     X_ji = gsl_matrix_get(X,j,i);
  //     X_ji = X_ji*X_ji;
  //     sum += X_ji;
  //   }
  //   sum = sqrt(sum);
  //   for(int j=0;j<n;j++)
  //   {
  //     X_ji = gsl_matrix_get(X,j,i);
  //     gsl_matrix_set(X,j,i,X_ji/sum);
  //   }
  // }
}

int main()
{
  gsl_matrix *Sm = gsl_matrix_alloc(n,n);
  gsl_matrix *Sdiag = gsl_matrix_alloc(n,n);
  gsl_matrix *Sevec = gsl_matrix_alloc(n,n);
  gsl_matrix *X = gsl_matrix_alloc(n,n);
  S(Sm,Sevec);
  matrixDiag(Sm,Sevec,Sdiag);
  diagonalizingMatrix(Sdiag,Sevec,X);

  gsl_matrix *tmp = gsl_matrix_alloc(n,n);
  gsl_matrix *isOne = gsl_matrix_alloc(n,n);
  gsl_matrix *X_T = gsl_matrix_alloc(n,n);
  gsl_matrix_transpose_memcpy(X_T,X);
  matrixProduct(Sm,X,tmp);
  printMatrix(20,tmp);
  printMatrix(4,X);
  matrixProduct(X_T,tmp,isOne);

  printMatrix(2,Sm);
  printMatrix(3,isOne);

  matrixProduct(X_T,X,tmp);
  printMatrix(20,tmp);
}
