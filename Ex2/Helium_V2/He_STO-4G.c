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
      gsl_matrix_set(Sm,p,q,tmp);
    }
  }
  return;
}

//Hamiltonian matrix (already integrated for GTO)
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
      A_pq = -2.*Z*pi/tmp; //Coulomb potential
      gsl_matrix_set(Hm,p,q,T_pq + A_pq);
    }
  }
  return;
}

//Direct matrix element POSSIBILE PROBLEMA: integrale inserito corretto?
double D_pq(int p, int q, int r, int s)
{
  double pi5 = pi*pi*pi*pi*pi;
  double tmp = alpha[p] + alpha[q] + alpha[r] + alpha[s];
  double V_pqrs = 2.*sqrt(pi5)/((alpha[p] + alpha[s])*(alpha[r] + alpha[q])*sqrt(tmp));

  return V_pqrs;
}

//Exchange matrix element POSSIBILE PROBLEMA: integrale inserito corretto?
double E_pq(int p, int q, int r, int s)
{
  double pi5 = pi*pi*pi*pi*pi;
  double tmp = alpha[p] + alpha[q] + alpha[r] + alpha[s];
  double V_pqrs = -2.*sqrt(pi5)/((alpha[p] + alpha[q])*(alpha[r] + alpha[s])*sqrt(tmp));

  return V_pqrs;
}

//Fock matrix
void F(gsl_matrix *Fm, gsl_matrix *Hm, gsl_matrix *Cm)
{
  int p,q,r,s,k;
  double H_pq, C_kr, C_ks, norm;
  double F_pq = 0;
  double P_rs = 0;
  norm = 0;

  for(p=0;p<n;p++)
  {
    for(q=0;q<n;q++)
    {
      F_pq=0;
      for(r=0;r<n;r++)
      {
        // if(r<p-1)
        // {
          for(s=0;s<n;s++)
          {
            P_rs = 0;
            for(k=0;k<=Z/2;k++)
            {
              C_kr = gsl_matrix_get(Cm,k,r);
              C_ks = gsl_matrix_get(Cm,k,s);
              P_rs = P_rs + C_kr*C_ks;//NOTE: sum over k up to Z/2
            }
            F_pq = F_pq + P_rs*(2*D_pq(p,q,r,s) + E_pq(p,q,r,s));
          }
        // }
        // else if(r==p-1)
        // {
          // for(s=0;s<q;s++)
          // {
          //   P_rs = 0;
          //   for(k=0;k<=Z/2;k++)
          //   {
          //     C_kr = gsl_matrix_get(Cm,k,r);
          //     C_ks = gsl_matrix_get(Cm,k,s);
          //     P_rs = P_rs + C_kr*C_ks;//sum over k up to N/2 (NOTA: quindi le ultime due colonne dei coefficienti non servono a niente?)
          //   }
          //   F_pq = F_pq + P_rs*(2*D_pq(p,q,r,s) + E_pq(p,q,r,s));
          // }
        // }
      }
      H_pq = gsl_matrix_get(Hm,p,q);
      F_pq = F_pq + H_pq;
      norm += F_pq*F_pq;
      // printf("p=%d, q=%d, F_pq=%f\n",p,q, F_pq);
      gsl_matrix_set(Fm,p,q,F_pq);
      gsl_matrix_set(Fm,q,p,F_pq);
    }
  }
  norm = 1./sqrt(norm);
  // gsl_matrix_scale(Fm,norm);
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
  tmp=0;
  double tmp2 = 0;
  for(p=0;p<n;p++)
  {
    tmp += gsl_vector_get(eval,p)*gsl_vector_get(eval,p);
    gsl_matrix_set(diagMat,p,p,gsl_vector_get(eval,p));
  }
  tmp = 1./sqrt(tmp);
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
  double X_ji,sum;
  for(int i=0;i<n;i++)
  {
    sum = 0;
    for(int j=0;j<n;j++)
    {
      X_ji = gsl_matrix_get(X,j,i);
      X_ji = X_ji*X_ji;
      sum += X_ji;
    }
    sum = sqrt(sum);
    for(int j=0;j<n;j++)
    {
      X_ji = gsl_matrix_get(X,j,i);
      gsl_matrix_set(X,j,i,X_ji/sum);
    }
  }
}

//F' = X_TFV, X = U_Ts^{-1/2}U, U = evec
void F_prime(gsl_matrix *Fm, gsl_matrix *Sm, gsl_matrix *evec, gsl_matrix *X)
{
  gsl_matrix *X_T = gsl_matrix_alloc(n,n);
  gsl_matrix *evec_T = gsl_matrix_alloc(n,n);
  gsl_matrix *Sdiag = gsl_matrix_alloc(n,n);
  gsl_matrix *tmp = gsl_matrix_alloc(n,n);

  matrixDiag(Sm,evec,Sdiag); // diagonalize Sm (resulting diagonal matrix is Sdiag)
  diagonalizingMatrix(Sdiag,evec,X); // build the matrix that diagonalizes Sm: X

  gsl_matrix_transpose_memcpy(X_T,X); // build the trnspose of X

  // build F'
  matrixProduct(Fm,X,tmp);
  matrixProduct(X_T,tmp,Fm);
  //normalize columns
  double X_ji,sum;
  for(int i=0;i<n;i++)
  {
    sum = 0;
    for(int j=0;j<n;j++)
    {
      X_ji = gsl_matrix_get(Fm,j,i);
      X_ji = X_ji*X_ji;
      sum += X_ji;
    }
    sum = sqrt(sum);
    for(int j=0;j<n;j++)
    {
      X_ji = gsl_matrix_get(Fm,j,i);
      gsl_matrix_set(Fm,j,i,X_ji/sum);
    }
  }
}

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

// Real symmetric eigensystem
void RSDE(gsl_matrix *Fm, gsl_matrix *Cm, gsl_vector *Ee)
{
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(Fm, Ee, Cm, w);
  gsl_eigen_symmv_free(w);
  return;
}

int main()
{
  gsl_matrix *Hm = gsl_matrix_alloc(n,n); // Hamiltonian matrix
  gsl_matrix *Sm = gsl_matrix_alloc(n,n); // Overlap matrix
  gsl_matrix *Sevec = gsl_matrix_alloc(n,n); // Overlap eigenvectors (from diagonalization of Sm)
  gsl_matrix *X = gsl_matrix_alloc(n,n); // Matrix that diagonalizes Sm
  gsl_matrix *Fm = gsl_matrix_alloc(n,n); // Fock matrix
  gsl_matrix *C_old = gsl_matrix_alloc(n,n); // Coefficient for basis expansion (each column is related to one eigenvalue)
  gsl_matrix *C_new = gsl_matrix_alloc(n,n); // Coefficients updated in the self-consistent procedure
  gsl_vector *Ee = gsl_vector_alloc(n); // Eigenvalues

  gsl_matrix_set_zero(Hm);
  gsl_matrix_set_zero(Sm);gsl_matrix_set_zero(Sevec);gsl_matrix_set_zero(X);
  gsl_matrix_set_zero(Fm);
  gsl_matrix_set_zero(C_old);
  gsl_matrix_set_zero(C_new);
  gsl_vector_set_zero(Ee);

  // First iteration (Fm = Hm)
  S(Sm,Sevec);
  H(Hm);
  F(Fm,Hm,C_new); //build Fock matrix
  printMatrix(0,Fm);
  F_prime(Fm,Sm,Sevec,X);
  printMatrix(1,Fm);
  printMatrix(2,Hm);
  RSDE(Fm,C_new,Ee); //Solve the generalized eigrnvalue problem
  printMatrix(3,C_new);
  gsl_matrix *tmp3 = gsl_matrix_alloc(n,n);
  matrixProduct(X,C_new,tmp3);
  gsl_matrix_memcpy(C_new,tmp3);

  double tmp = fabs(minE(Ee));//POSSIBLE PROBLEM IF ALREADY SMALLER THAN convergence
  printf("minE = %f\n", tmp);

  //self-consistent procedure
  double alfa = 5e-1; //coefficient for the mixture of the old and new coefficients ([10^-3, 0.5])
  printf("alfa = %e\n", alfa);
  double check = 1; // store the variation of energy (minimum) between cycles
  double convergence = 1e-10; // cycle stops when the energy difference (check) is lower than "convergence"
  double evec[n][n], eval[n]; //store eigenvectors and eigenvalues

  double tmp2,sum,C_kl,H_kl,F_kl;
  int i,j,k;
  k=0;
  while(check>convergence)
  {
    //compute new coefficients as a small variation of old coefficients
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        tmp2 = alfa*gsl_matrix_get(C_new,i,j) + (1-alfa)*gsl_matrix_get(C_old,i,j);
        gsl_matrix_set(C_new,i,j,tmp2);
      }
    }
    // Recalculate S and H matrix (RSDE modifies them)
    F(Fm,Hm,C_new); //build Fock matrix
    F_prime(Fm,Sm,Sevec,X);
    gsl_matrix_memcpy(C_old,C_new);//update C_old
    RSDE(Fm,C_new,Ee); //Solve the generalized eigrnvalue problem
    matrixProduct(X,C_new,tmp3);
    gsl_matrix_memcpy(C_new,tmp3);
    // printMatrix(k+4,C_new);
    sum = 0;
    // calculate energy
    for(int k=0;k<n;k++)
    {
      for(int l=0;l<n;l++)
      {
        C_kl = gsl_matrix_get(C_new,k,l);
        H_kl = gsl_matrix_get(Hm,k,l);
        F_kl = gsl_matrix_get(Fm,k,l);
        sum = C_kl*(H_kl + F_kl);
      }
    }
    tmp2 = fabs(minE(Ee));
    check = fabs(tmp2-tmp);
    tmp = tmp2;
    printf("----------------------%d-------------------------\n", k);
    printf("minE = %f       check = %e\n", minE(Ee), check);
    printf("minE = %f       check = %e\n", sum, check);
    k++;
    // if(k>150) check = 0;
  }

  printf("------ minE = %f ------\n", minE(Ee));

  return 0;
}
