
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

    if(posNew[i].x >= 2*a) posNew[i].x -= 2*a;
    if(posNew[i].x < 0)    posNew[i].x += 2*a;

    if(posNew[i].y >= 2*a) posNew[i].y -= 2*a;
    if(posNew[i].y < 0)    posNew[i].y += 2*a;

    if(posNew[i].z >= 2*a) posNew[i].z -= 2*a;
    if(posNew[i].z < 0)    posNew[i].z += 2*a;
  }
}

double r_ij(struct point points[N], int i, int j)
{
  double dx = fabs(points[i].x - points[j].x);
  double dy = fabs(points[i].y - points[j].y);
  double dz = fabs(points[i].z - points[j].z);

  if(dx > a) dx = 2*a - dx;
  if(dy > a) dy = 2*a - dy;
  if(dz > a) dz = 2*a - dz;

  return gsl_hypot3(dx, dy, dz);
}

int MC_acc()
{
  double acc;

  double sumPsi = 0;
  for (int i = 0; i < N; i++)
  for (int j = i+1; j < N; j++)
  {
    sumPsi += gsl_pow_5(1/r_ij(pos, i, j)) - gsl_pow_5(1/r_ij(posNew, i, j));
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

double normScalarProd(struct point points[N], int l, int i, int j)
{
  double prod = (points[l].x - points[i].x) * (points[l].x - points[j].x) +
                (points[l].y - points[i].y) * (points[l].y - points[j].y) +
                (points[l].z - points[i].z) * (points[l].z - points[j].z);

  return prod/(r_ij(posNew, l, i) * r_ij(posNew, l, j));
}

double E_kin()
{
  double Etmp = 0;
  double Etmp2 = 0;

  for (int l = 0; l < N; l++)
  {
    for (int i = 0; i < N; i++)
    {
      if(l != i)
      {
        Etmp += gsl_pow_7(1/r_ij(pos, i, l));

        for (int j = 0; j < N; j++)
        {
          if(l != j)
          {
            Etmp2 += normScalarProd(pos, l, i, j) * gsl_pow_6(1/r_ij(pos, l, i)) * gsl_pow_6(1/r_ij(pos, l, j));
          }
        }
      }
    }
  }
  return h2_2m * (25/4* b10 * Etmp2 + 10 * b5 * Etmp);
}

double E_kin2()
{
  double Etmp = 0;
  double Etmp2 = 0;

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    if(i != j)
    {
      Etmp += gsl_pow_7(1/r_ij(pos, i, j));
    }
  }
  return h2_2m * 5 * b5 * Etmp;
}

double E_pot()
{
  double Etmp = 0;

  for (int i = 0; i < N; i++)
  for (int j = i+1; j < N; j++)
  {
    Etmp += 4*(gsl_pow_int(r_ij(pos, i, j), -12) - gsl_pow_int(r_ij(pos, i, j), -6));
  }

  return 0;
}


int main()
{
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  a = 1;
  Delta = a* 0.007;
  b = 1;
  b5 = pow(b, 5);
  b10 = b5*b5;

  placeParticles();

  FILE *fp;
  fp=fopen("ThermAcc.csv","w+");

  for (int t = 0; t < 1000; t++)
  {

    int tot = 0;
    for (int tt = 0; tt < 200; tt++)
    {
      MC_Move();
      tot += MC_acc();

      //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
      if(t > 500)
      {
          Nsamples++;
          Eloc += E_kin(); // + E_pot(); //calculate the energy
          Eloc2 += E_kin2(); // + E_pot(); //calculate the energy
      }
    }

    //print the progress bar
    if(t % 5 == 0) printf("%d\n", t/5);

    //print data into the file
    fprintf(fp,"%d %f\n", t, tot/200.0);
  }

  fclose(fp);

  printf("%lg \t %lg\n", Eloc/Nsamples, Eloc2/Nsamples);

  writeCSVparticles(filename, pos);
  writeCSVparticles(filename1, posNew);

  return 0;
}
