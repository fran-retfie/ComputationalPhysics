#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

int n=3;
void printMatrix(char matrixName, gsl_matrix *m)
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
  printf("----------------------------------\n");
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

int main()
{
  //matrix to diagonalize
  gsl_matrix *A = gsl_matrix_alloc(3,3);
  gsl_matrix_set(A,0,0,1);gsl_matrix_set(A,0,1,1);gsl_matrix_set(A,0,2,3);
  gsl_matrix_set(A,1,0,1);gsl_matrix_set(A,1,1,5);gsl_matrix_set(A,1,2,1);
  gsl_matrix_set(A,2,0,3);gsl_matrix_set(A,2,1,1);gsl_matrix_set(A,2,2,1);

  //identity
  gsl_matrix *I = gsl_matrix_alloc(3,3);
  gsl_matrix_set_identity(I);

  char a = {'A'};
  printMatrix(a,A);

  matrixProduct(A,A,I);
  char i = {'I'};
  printMatrix(i,I);

  //eigenvalues and eigenvectors
  gsl_vector *eval = gsl_vector_alloc(3);
  gsl_matrix *evec = gsl_matrix_alloc(3,3);
  gsl_matrix *evecT = gsl_matrix_alloc(3,3);
  gsl_matrix *D = gsl_matrix_alloc(3,3);
  gsl_vector_set_zero(eval);
  gsl_matrix_set_zero(evec);
  gsl_matrix_set_zero(evecT);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(A,eval,evec,w);
  gsl_matrix_transpose_memcpy(evecT, evec);
  printMatrix(a,evecT);
  gsl_matrix_set(D,0,0,gsl_vector_get(eval,0));gsl_matrix_set(D,0,1,0);gsl_matrix_set(D,0,2,0);
  gsl_matrix_set(D,1,0,0);gsl_matrix_set(D,1,1,gsl_vector_get(eval,1));gsl_matrix_set(D,1,2,0);
  gsl_matrix_set(D,2,0,0);gsl_matrix_set(D,2,1,0);gsl_matrix_set(D,2,2,gsl_vector_get(eval,2));
  printMatrix(a,D);
  matrixProduct(D,evec,D);
  matrixProduct(evecT,D,D);


  printMatrix(a,D);

  return 0;
}
