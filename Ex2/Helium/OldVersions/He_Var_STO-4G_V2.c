#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

double pi=4.*atan(1.);
int n=4;//# of basis functions
int Z=2;//# of electrons
//gaussians exponents
double alpha[4] = {14.899983, 2.726485, 0.757447, 0.251390};

//gaussian basis (unused)
double Chi(double r, double alpha)
{
  double r2 = r*r;

  return (exp(-alpha*r2));
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

  //diagonalization
  // gsl_vector *eval = gsl_vector_alloc(n);
  // gsl_vector_set_zero(eval);
  // gsl_matrix_set_zero(evec);
  // gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  // gsl_eigen_symmv(Sm,eval,evec,w);
  // gsl_matrix_set_zero(Sm);
  // tmp=1;
  // for(p=0;p<n;p++)
  // {
  //   tmp = tmp*gsl_vector_get(eval,p);
  //   gsl_matrix_set(Sm,p,p,gsl_vector_get(eval,p));
  // }
  //
  // tmp = 1./sqrt(tmp);
  //
  // //normalization
  // gsl_matrix_scale(Sm,tmp);

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

// double doubleInt()


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
      F_pq = 0;
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
            printf("%f\n", P_rs);
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

//F' = V_TFV, V = S^{-1/2}U, U = Sevec
void F_prime(gsl_matrix *Fm, gsl_matrix *Sm, gsl_matrix *Sevec, gsl_matrix *V)
{
  gsl_matrix *V_T = gsl_matrix_alloc(n,n);
  gsl_matrix_memcpy(V,Sevec);

  //build V matrix
  int i,j;
  double lambda_i,tmp;
  for(i=0;i<n;i++)
  {
    lambda_i = sqrt(gsl_matrix_get(Sm,i,i));
    tmp = lambda_i*gsl_matrix_get(V,i,i);
    gsl_matrix_set(V,i,i,tmp);
  }

  gsl_matrix_transpose_memcpy(V_T,V);//V_T = transpose of V

  matrixProduct(Fm,V,Fm);//Fm updated to be FV
  matrixProduct(V_T,Fm,Fm);//Fm updated to V_TFV = F'
}

// find minimum component of a vector (used for energy minimum)
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
  gsl_matrix *Sevec = gsl_matrix_alloc(n,n);
  gsl_matrix *Fm = gsl_matrix_alloc(n,n);
  gsl_matrix *Vm = gsl_matrix_alloc(n,n);
  gsl_matrix *C_old = gsl_matrix_alloc(n,n);
  gsl_matrix *C_new = gsl_matrix_alloc(n,n);
  gsl_vector *Ee = gsl_vector_alloc(n);
  gsl_matrix_set_zero(Hm);
  gsl_matrix_set_zero(Sm);
  gsl_matrix_set_zero(Fm);
  gsl_matrix_set_zero(C_old);
  gsl_matrix_set_zero(C_new);
  gsl_vector_set_zero(Ee);

  //self-consistent procedure
  double alfa = 1e-2;//coefficient in the variation of the coefficients
  double tmp,norm,tmp2;
  int i,j;

  //first iteration (to set the min of energy)
  S(Sm, Sevec); //build overlap matrix
  char s={'S'};char h={'H'};char c_n={'C'};char c_o={'O'};char f={'F'};
  printMatrix(s,Sm);
  H(Hm); //build Hamiltonian
  F(Fm,Hm,C_new); //build Fock matrix
  // F_prime(Fm,Sm,Sevec,Vm);
  RGSDE(Fm,Sm,C_new,Ee); //Solve the generalized eigrnvalue problem
  tmp = fabs(minE(Ee));//POSSIBLE PROBLEM IF ALREADY SMALLER THAN convergence
  double check = 1;
  double convergence = 1e-6;
  printf("minE = %f       check = %f\n", minE(Ee), check);
  while(check>convergence)
  {
    //compute new coefficients as a small variation of
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
      {
        tmp2 = alfa*gsl_matrix_get(C_new,i,j) + (1-alfa)*gsl_matrix_get(C_old,i,j);
        // printf("tmp2 = %f\n", tmp2);
        gsl_matrix_set(C_new,i,j,tmp2);
        //norm = gsl_matrix_norm1(C_new);
      }
    }
    S(Sm,Sevec); //build overlap matrix
    H(Hm); //build Hamiltonian

    gsl_matrix_memcpy(C_old,C_new);//update C_old
    F(Fm,Hm,C_new); //build Fock matrix
    // F_prime(Fm,Sm,Sevec,Vm);
    // printMatrix(f,Fm);
    RGSDE(Fm,Sm,C_new,Ee); //Solve the generalized eigrnvalue problem
    tmp2 = fabs(minE(Ee));
    check = fabs(tmp2-tmp);
    tmp = tmp2;
    // check = 0;
    printf("minE = %f       check = %f\n", minE(Ee), check);
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
