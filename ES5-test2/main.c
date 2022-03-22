
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
double b10;
double Eloc = 0;
double Eloc2 = 0;
double Var = 0;
double Var2 = 0;

long Nsamples = 0;
long tot = 0;

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

void MC_Move()
{
  for(int i=0;i<N;i++)
  {
    posNew[i].x = pos[i].x + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].y = pos[i].y + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].z = pos[i].z + gsl_ran_flat(r, -Delta, Delta);
  }
}

double r_ij(struct point points[N], int i)
{
  double dx = points[i].x;
  double dy = points[i].y;
  double dz = points[i].z;

  return gsl_hypot3(dx, dy, dz);
}

int MC_acc()
{
  double acc;

  double sumPsi = 0;
  for (int i = 0; i < N; i++)
  {
    sumPsi += r_ij(pos, i) - r_ij(posNew, i);
  }

  acc = exp(2*b*sumPsi);

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

double E_L()
{
  double E = 0;
  for (int i = 0; i < N; i++) E += -1/r_ij(pos, i) - b/2*(b-2/r_ij(pos, i));
  return E;
}

double E_L2()
{
  double E = 0;
  for (int i = 0; i < N; i++) E += -1/r_ij(pos, i) + b/4/r_ij(pos, i);
  return E;
}

int main()
{
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  a = 1;
  Delta = a * 0.1;
  b = 1;
  b5 = pow(b, 5);
  b10 = b5*b5;

  double eloc;
  double eloc2;
  double var;
  double var2;

  placeParticles();

  FILE *fp;
  fp=fopen("ThermAcc.csv","w+");

  for (int t = 0; t < 100000; t++)
  {

    MC_Move();
    tot += MC_acc();

    //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
    if(t > 500)
    {
        Nsamples++;
        Eloc += E_L(); //calculate the energy
        Var += E_L()*E_L(); //calculate the energy
        Eloc2 += E_L2();
        Var2 += E_L2()*E_L2();

    }

    //print the progress bar
    //if(t % 5 == 0) printf("%d\n", t/5);

    //print data into the file
    fprintf(fp,"%d %f\n", t, tot/200.0);
  }

  fclose(fp);

  eloc = Eloc/Nsamples;
  eloc2 = Eloc2/Nsamples;
  var = Var/Nsamples;
  var2 = Var2/Nsamples;
  printf("%lg \t %lg \n", eloc, sqrt((var - eloc*eloc)));
  printf("%lg \t %lg \n", eloc2, sqrt((var2 - eloc2*eloc2)));

  printf("%lg\n", (double) tot/(double) Nsamples);

  return 0;
}
