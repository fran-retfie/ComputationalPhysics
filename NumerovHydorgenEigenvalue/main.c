
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

#include "CSVio.h"
#include "DiffInt.h"

#ifndef long_double_t
  typedef  long double long_double_t; //__float128  long_double_t;
#endif

#define N 1000
#define Nmax 10
#define Lmax 4

const double xStart = 0.01;
const double xEnd = 20;
const double h = (xEnd-xStart)/(double) N;

const double C = 13.605693122904;

double Ea[2][Nmax];
double psi[2][N];

int n = 1;
int l = 0;
double E = 0;

//effective potential (normalized!!)
double fr(double x)
{
  if(x == 0) printf("error\n");
  double r1 = 1/x;
  double pot = -r1;
  if(l != 0) pot += l*(l+1)*r1*r1*0.5;
  return pot;
}

long_double_t CalcNorm()
{
  long_double_t I = 0;
  for(int j = 1; j < N-1; j+=2)
  {
    I += h/3*(psi[1][j-1]*psi[1][j-1]+4*psi[1][j]*psi[1][j]+psi[1][j+1]*psi[1][j+1]);
  }
  return sqrt(I);
}

gsl_matrix *invert_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(N);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(N,N);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

int main()
{

  for(l = 0; l<Lmax; l++)
  {
    printf("\n\nl = %i\n", l);

    const double inv12 = 1/12.0;
    double h2 = 1/(h*h);
    double x;

    gsl_vector *xvec = gsl_vector_alloc(N);
    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_matrix *evec = gsl_matrix_alloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_matrix *Bi = gsl_matrix_calloc(N, N);
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *H = gsl_matrix_calloc(N, N);

    for(int i=0;i<N;i++)
    {
      //write x vector
      x = xStart+i*h;
      gsl_vector_set(xvec,i,x);

      //write B matrix
      gsl_matrix_set(B, i, i, inv12*10);
      if(i<N-1) gsl_matrix_set(B, i, i+1, inv12);
      if(i>0)   gsl_matrix_set(B, i, i-1, inv12);

      //write A matrix
      if((i!=0) && (i!=(N-1)))  gsl_matrix_set(A, i, i, -2*h2);
      if(i<N-1) gsl_matrix_set(A, i, i+1, h2);
      if(i>0)   gsl_matrix_set(A, i, i-1, h2);

      //write V matrix
      gsl_matrix_set(H, i, i, fr(x));
    }

    //calculate H = B^-1 A + V
    Bi = invert_matrix(B);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -0.5, A, Bi, 1.0, H);

/*
    for(int a = 0; a<N; a++)
    {
      for(int b = 0; b<N; b++) printf("%2.2g\t", gsl_matrix_get(H, a, b));
      printf("\n");
    }
    printf("\n\n\n\n");
*/

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);
    gsl_eigen_symmv(H, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

/*
    for(int a = 0; a<N; a++)
    {
      for(int b = 0; b<N; b++) printf("%2.3g\t", gsl_matrix_get(evec, a, b));
      printf("\n");
    }
*/

    int i = 0;
    for(int n = 0; n < Nmax; i++)
    {
      for(int j = 0; j < N; j++)
      {
        psi[0][j] = gsl_vector_get(xvec, j);
        psi[1][j] = gsl_matrix_get(evec, j, n);
      }

      double E = gsl_vector_get(eval, i)*2*C;
      if(E>-20)
      {
        printf("n = %i, E = %15.14g, Err = %15.14g\n",n,E, (E+C)/C);
        Ea[1][n] = E;

        double norm = CalcNorm();
        for(int j = 0; j < N; j++) psi[1][j] = psi[1][j]/norm;

        char filename[13];
        sprintf(filename,"dati/plot%02i%02i",n,l);
        writeCSVdouble(filename,(double *)psi, 2, N);

        n++;
      }
    }
  }

  //gsl_vector_free(eval);
  //gsl_matrix_free(evec);

  writeCSVdouble("energies.csv",(double *)Ea, 2, Nmax);

  return 0;
}
