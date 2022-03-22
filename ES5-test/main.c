
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

long Nsamples = 0;
long tot = 0;
struct point pos;
struct point posNew;

void placeParticles()
{
  pos.x = 0.5;
  pos.y = 1.0;
  pos.z = 2.0;
}

void MC_Move()
{
  posNew.x = pos.x + gsl_ran_flat(r, -Delta, Delta);
  posNew.y = pos.y + gsl_ran_flat(r, -Delta, Delta);
  posNew.z = pos.z + gsl_ran_flat(r, -Delta, Delta);
}

double r_ij(struct point points)
{
  return gsl_hypot3(points.x, points.y, points.z);
}

int MC_acc()
{

  double acc = exp(-2*b*(r_ij(posNew) - r_ij(pos))); //(r_ij(posNew)/r_ij(pos))*(r_ij(posNew)/r_ij(pos)) *

  if(acc > gsl_ran_flat(r, 0, 1))
  {
    pos.x = posNew.x;
    pos.y = posNew.y;
    pos.z = posNew.z;
    return 1;
  }
  else  return 0;
}


double E_L()
{
  return -1/r_ij(pos) - b/2*(b-2/r_ij(pos));
}

double eloc;
double eloc2;

int main()
{
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  a = 1;
  Delta = a * 2;
  b = 1.1;

  placeParticles();

  long tot = 0;
  for (int t = 0; t < 100000; t++)
  {
    MC_Move();
    tot += MC_acc();

    //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
    if(t > 500)
    {
        Nsamples++;
        printf("%lg \t %lg \n", r_ij(pos), r_ij(posNew));
        Eloc += E_L(); // + E_pot(); //calculate the energy
        Eloc2 += E_L()*E_L(); // E_kin2(); // + E_pot(); //calculate the energy
    }

    //print data into the file
    //fprintf(fp,"%d %f\n", t, tot/200.0);
  }

  eloc = Eloc/Nsamples;
  eloc2 = Eloc2/Nsamples;
  printf("%lg \t %lg \n", eloc, sqrt((eloc2 - eloc*eloc)));

  printf("%lg\n", (double) tot/(double) Nsamples);

  return 0;
}
