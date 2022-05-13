#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

double pi=4.*atan(1.);
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
void FockMatrix(gsl_matrix *Fm, double Cm[N][n], double Dm[n][n][n][n], double Hm[n][n])
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
            F_pq += Cm[k][r]*Cm[k][s]*(2*Dm[p][r][q][s] - Dm[p][r][s][q]);
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
void E_HF(double Cm[N][n], double eval[orbitalNum], double Dm[n][n][n][n], double Hm[n][n], double E[2])
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
              tmp += Cm[i][p]*Cm[k][r]*Cm[i][q]*Cm[k][s]*(2*Dm[p][r][q][s]-Dm[p][r][s][q]); //No 2 for H
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
void updateCoef(double Cm[N][n], double CmNew[N][n])
{
  for(int k=0;k<orbitalNum;k++)//column
  {
    for(int i=0;i<n;i++)//row
    {
      Cm[k][i] = mixing*CmNew[k][i] + (1-mixing)*Cm[k][i];
    }
  }
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
  //Matrices init
  double Hm[n][n];
  gsl_matrix *Sm = gsl_matrix_alloc(n,n);
  double Dm[n][n][n][n];
  double Cm[N][n];
  double CmNew[N][n];
  double Eval[orbitalNum];

  gsl_matrix *Fm = gsl_matrix_alloc(n,n);
  gsl_matrix *evec = gsl_matrix_alloc(n,n);
  gsl_vector *eval = gsl_vector_alloc(n);
  gsl_matrix_set_zero(evec);
  gsl_matrix_set_zero(Fm);
  gsl_vector_set_zero(eval);

  double thr = 0.00001;
  int counter = 0;
  double E[2];
  E[0]=0;E[1]=1;

  H(Hm);
  S(Sm);
  D(Dm);
  initCoef(Cm);
  printCMatrix(Cm);

  printDoubleMatrix(h,Hm);
  printGslMatrix(s,Sm);
  //printDMatrix(de,Dm);

  double tmp,tmpOld,Eold;
  int index[orbitalNum];
  Eold=1;

  while(fabs(E[0]-Eold)>thr)
  //for(int ii=0;ii<200;ii++)
  {
    counter++;
    for(int k=0;k<orbitalNum;k++)
    {
      index[k]=0;
      Eval[k]=0;
    }
    S(Sm);
    tmp = 0;
    tmpOld = 0;
    FockMatrix(Fm,Cm,Dm,Hm);
    gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(Fm, Sm, eval, evec, w);
    gsl_matrix_set_zero(Fm);
    //printGslMatrix(h,evec);
    //Sort eval and evec
    for(int k=0;k<orbitalNum;k++)
    {
      for(int i=0;i<n;i++)
      {
        tmp = gsl_vector_get(eval,i);
        if((tmp<tmpOld) && ((i!=index[k-1])||(k==0)))
        {
          Eval[k] = tmp;
          index[k] = i;
          tmpOld = tmp;
        }
      }
    }
    //Get new coeff from evec
    for(int k=0;k<orbitalNum;k++)//column
    {
      for(int i=0;i<n;i++)//row
      {
        tmp = gsl_matrix_get(evec,i,index[k]);
        CmNew[k][i] = tmp;
      }
    }

    updateCoef(Cm,CmNew);
    Eold = E[0];
    E_HF(Cm,Eval,Dm,Hm,E);
    printf("E1: %lf     E2: %lf\n", E[0],E[1]);

    gsl_eigen_gensymmv_free(w);
  }
  printf("E: %lf    iterations: %d\n", E[0], counter);
  printCMatrix(Cm);
  //printCMatrix(CmNew);
  tmp=0;
  for(int i=0;i<n;i++)
  {
    tmp += Cm[0][i]*Cm[0][i];
  }
  printf("%lf\n", sqrt(tmp));
  printf("Cm matrix: \n");
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<n;j++)
    {
      Cm[i][j] = Cm[i][j]/sqrt(tmp);
      printf("%f  ", Cm[i][j]);
    }
    printf("\n");
  }
  tmp=0;
  for(int i=0;i<n;i++)
  {
    tmp += Cm[0][i]*Cm[0][i];
  }
  printf("%lf\n", sqrt(tmp));
}
