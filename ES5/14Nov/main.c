
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "CSVio.h"

#define N 32

/* Declare the struct with coordinates of a point x, y, z */
struct point {
   double    x;
   double    y;
   double    z;
};

const gsl_rng_type * T;
gsl_rng * r;

char filename[20] = "initialPos.csv";
char filename1[20] = "newPos.csv";

int const Ncells = 2;
double a;
double Delta;
double b;
double b5;
double Eloc = 0;
double Eloc2 = 0;

long Nsamples = 0;

double const h2_2m = 1;

// Position of all particles
struct point pos[N];
struct point posNew[N];

void placeParticles()
{
  int d = 0;

  for (int i = 0; i < Ncells; i++)
  for (int j = 0; j < Ncells; j++)
  for (int k = 0; k < Ncells; k++)
  {
    pos[d].x = a*(i + 0);
    pos[d].y = a*(j + 0);
    pos[d].z = a*(k + 0);

    pos[d+1].x = a*(i + 0.5);
    pos[d+1].y = a*(j + 0.5);
    pos[d+1].z = a*(k + 0);

    pos[d+2].x = a*(i + 0);
    pos[d+2].y = a*(j + 0.5);
    pos[d+2].z = a*(k + 0.5);

    pos[d+3].x = a*(i + 0.5);
    pos[d+3].y = a*(j + 0);
    pos[d+3].z = a*(k + 0.5);

    d += 4;
  }
}

int writeCSVparticles(char *filename, struct point points[N])
{
  FILE *fp;
  fp=fopen(filename,"w+");

  for(int i=0;i<N;i++)
  {
      fprintf(fp,"%lg ", points[i].x );
      fprintf(fp,",%lg ", points[i].y );
      fprintf(fp,",%lg ", points[i].z );
      fprintf(fp,",%i \n", i );
  }

  fclose(fp);
  return 0;
}

struct point PBC(struct point dr)
{
  dr.x = dr.x - 2*a*rint(dr.x/(a));
  dr.y = dr.y - 2*a*rint(dr.y/(a));
  dr.z = dr.z - 2*a*rint(dr.z/(a));

  return dr;
}

void MC_Move()
{
  for(int i=0;i<N;i++)
  {
    posNew[i].x = pos[i].x + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].y = pos[i].y + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].z = pos[i].z + gsl_ran_flat(r, -Delta, Delta);

    if(posNew[i].x >= 2*a) posNew[i].x -= 2*a;
    if(posNew[i].x < 0)    posNew[i].x += 2*a;

    if(posNew[i].y >= 2*a) posNew[i].y -= 2*a;
    if(posNew[i].y < 0)    posNew[i].y += 2*a;

    if(posNew[i].z >= 2*a) posNew[i].z -= 2*a;
    if(posNew[i].z < 0)    posNew[i].z += 2*a;
  }
}

double r_il(struct point dr)
{
  return gsl_hypot3(dr.x, dr.y, dr.z);
}

int MC_acc()
{
  double acc;
  struct point dr;
  double sumPsi = 0;

  for (int i = 0; i < N; i++)
  {
    for (int j = i+1; j < N; j++)
    {
      dr.x = pos[i].x - pos[j].x;
      dr.y = pos[i].y - pos[j].y;
      dr.z = pos[i].z - pos[j].z;
      dr = PBC(dr);
      sumPsi += gsl_pow_int(r_il(dr), -5) - gsl_pow_int(r_il(dr), -5);
    }
  }
  acc = exp(b5*sumPsi);

  if(acc > gsl_ran_flat(r, 0, 1))
  {
    for (int i = 0; i < N; i++)
    {
      pos[i].x = posNew[i].x;
      pos[i].y = posNew[i].y;
      pos[i].z = posNew[i].z;
    }
    return 1;
  }
  else return 0;

}

double E_kin()
{
  double Ek = 0;
  struct point dr;
  double der = 0;
  double du1_l = 0;
  struct point du2[N];

  for (int l = 0; l < N; l++)
  {
    du1_l = 0;
    du2[l].x = 0;du2[l].y = 0;du2[l].z = 0;
    //printf("******************************************\n");
    for (int i = 0; i < N; i++)
    {
      if(i != l)
      {
        //printf("%d %d\n", l,i);
        //printf("pre l = %d   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",l, dr.x,dr.y,dr.z);
        dr.x = pos[l].x - pos[i].x;
        dr.y = pos[l].y - pos[i].y;
        dr.z = pos[l].z - pos[i].z;
        dr = PBC(dr);
        //printf("post l = %d   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",l, dr.x,dr.y,dr.z);
        der = 5*b5*gsl_pow_int(r_il(dr), -7); //derivative times 1/r_il
        du1_l += 2*der;
        du2[l].x += der*dr.x;
        du2[l].y += der*dr.y;
        du2[l].z += der*dr.z;
        //printf("der: %lg\n", der);
        //printf("l = %d   du2[l].x = %lg   du2[l].y = %lg    du2[l].z = %f \n",l, du2[l].x,du2[l].y,du2[l].z);
      }
      //printf("l = %d   du2[l].x = %lg   du2[l].y = %lg    du2[l].z = %lg \n",l, du2[l].x,du2[l].y,du2[l].z);
    }
    Ek += du1_l - du2[l].x*du2[l].x - du2[l].y*du2[l].y - du2[l].z*du2[l].z;
    //printf("l = %d   Ek = %lg\n",l, Ek);
  }
  return -h2_2m * Ek;
}

double E_kinJF()
{
  double Etmp = 0;
  struct point dr;

  for (int l = 0; l < N; l++)
  {
    for (int i = 0; i < N; i++)
    {
      if(i != l)
      {
        dr.x = pos[l].x - pos[i].x;
        dr.y = pos[l].y - pos[i].y;
        dr.z = pos[l].z - pos[i].z;
        dr = PBC(dr);
        Etmp += 5*b5*gsl_pow_int(r_il(dr), -7);

      }
    }
  }
  return h2_2m * Etmp;
}

double E_pot()
{
  double Etmp = 0;
  struct point dr;

  for (int i = 0; i < N; i++)
  {
    for (int j = i+1; j < N; j++)
    {
      dr.x = pos[i].x - pos[j].x;
      dr.y = pos[i].y - pos[j].y;
      dr.z = pos[i].z - pos[j].z;
      dr = PBC(dr);
      Etmp += 4*(gsl_pow_int(r_il(dr), -12) - gsl_pow_int(r_il(dr), -6));
    }
  }
  return Etmp;
}

void test(int i, int l)
{
  struct point dr;
  double der = 0;

  dr.x = pos[l].x - pos[i].x;
  dr.y = pos[l].y - pos[i].y;
  dr.z = pos[l].z - pos[i].z;
  printf("prev   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",dr.x,dr.y,dr.z);
  dr = PBC(dr);
  printf("post   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",dr.x,dr.y,dr.z);
  der = 5*b5*gsl_pow_int(r_il(dr), -7); //derivative times 1/r_il

  printf("der: %lg\n", der);
}

int main()
{
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  a = 5.;
  Delta = a* 0.002;
  b = 1.;
  b5 = pow(b, 5);

  placeParticles();
/*
  double der = 0;
  struct point dr;
  struct point du2[N];
  for(int l=0; l<N;l++)
  {
    for(int i=0; i<N;i++)
    {
      if(i != l)
      {
        printf("l:%d\n", l);
        test(i,l);
        dr.x = pos[l].x - pos[i].x;
        dr.y = pos[l].y - pos[i].y;
        dr.z = pos[l].z - pos[i].z;
        printf("prev   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",dr.x,dr.y,dr.z);
        dr = PBC(dr);
        printf("post   dr.x = %lg   dr.y = %lg    dr.z = %lg \n",dr.x,dr.y,dr.z);
        der = 5*b5*gsl_pow_int(r_il(dr), -7);
        du2[l].x += der*der*dr.x*dr.x;
        du2[l].y += der*der*dr.y*dr.y;
        du2[l].z += der*der*dr.z*dr.z;
        printf("der: %lg\n", der);
        printf("l = %d   du2[l].x = %lg   du2[l].y = %lg    du2[l].z = %f \n",l, du2[l].x,du2[l].y,du2[l].z);
      }
    }
  }
*/

  double E_k = E_kin();
  double E_p = E_pot();
  double E_kJF = E_kinJF();
  printf("E_kin: %lg    E_kinJF: %lg     E_pot: %lg\n", E_k, E_kJF, E_p);


  FILE *fp;
  fp=fopen("ThermAcc.csv","w+");

  for (int t = 0; t < 500; t++)
  {

    int tot = 0;
    for (int tt = 0; tt < 1000; tt++)
    {
      MC_Move();
      tot += MC_acc();

      //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
      if(t > 150)
      {
          Nsamples++;
          E_p += E_pot();
          Eloc += E_kin();
          Eloc2 += E_kinJF();
      }
    }

    fprintf(fp,"%d %d\n", t, tot);
  }

  fclose(fp);

  printf("%lg \t %lg\n", Eloc/Nsamples, Eloc2/Nsamples);

  writeCSVparticles(filename, pos);
  writeCSVparticles(filename1, posNew);

  return 0;
}
