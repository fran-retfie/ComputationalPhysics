
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <time.h>

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

int const Ncells = 2;
double a;
double Delta;
double b;
double b5;
double b10;
double Eloc;
double Eloc2;
double Var;
double Var2;
double pot;
double pot2;
double kin;
double kin2;
double P1 = 1;
double P2 = 1;

long Nsamples;
long tot;

struct point point_allias;

double const h2_2m = 0.09075588736790198;

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

void compute_allias(point& me, point& oth)
{
  if((me.x - oth.x) >  a) oth.x += 2*a;
  if((me.y - oth.y) >  a) oth.y += 2*a;
  if((me.z - oth.z) >  a) oth.z += 2*a;
  if((me.x - oth.x) < -a) oth.x -= 2*a;
  if((me.y - oth.y) < -a) oth.y -= 2*a;
  if((me.z - oth.z) < -a) oth.z -= 2*a;
}

double ran_flat(double Delta)
{
  return Delta * ( 2*(double) rand()/(double) RAND_MAX - 1 );
}

void MC_Move()
{
  for(int i=0;i<N;i++)
  {
    posNew[i].x = pos[i].x + ran_flat(Delta);
    posNew[i].y = pos[i].y + ran_flat(Delta);
    posNew[i].z = pos[i].z + ran_flat(Delta);

    if(posNew[i].x > a)  posNew[i].x -= 2*a;
    if(posNew[i].y > a)  posNew[i].y -= 2*a;
    if(posNew[i].z > a)  posNew[i].z -= 2*a;
    if(posNew[i].x < -a) posNew[i].x += 2*a;
    if(posNew[i].y < -a) posNew[i].y += 2*a;
    if(posNew[i].z < -a) posNew[i].z += 2*a;
    //posNew[i].x += 2*a*rint(posNew[i].x/(2*a));
    //posNew[i].y += 2*a*rint(posNew[i].y/(2*a));
    //posNew[i].z += 2*a*rint(posNew[i].z/(2*a));

  }
}

double r_ij(struct point points[N], point& point_i, int j)
{
  double dx = (points[j].x - point_i.x);
  double dy = (points[j].y - point_i.y);
  double dz = (points[j].z - point_i.z);

  return gsl_hypot3(dx, dy, dz);
}

int MC_acc()
{
  double sumPsi1 = 0;
  double sumPsi2 = 0;

  for (int i = 0;   i < N; i++)
  for (int j = i+1; j < N; j++)
  {
    point_allias = posNew[i];
    compute_allias(posNew[j], point_allias);
    sumPsi1 += pow(1/r_ij(posNew, point_allias, j), 5);

    point_allias = pos[i];
    compute_allias(pos[j], point_allias);
    sumPsi2 += pow(1/r_ij(pos, point_allias, j), 5);
  }

  double acc = exp(b5*(sumPsi2 - sumPsi1));

  double randN = (double) rand()/(double) RAND_MAX; //gsl_ran_flat(r, 0, 1);
  if(acc > randN)
  {
    //printf("%lg %lg\n",acc ,randN);
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
  double Etmp = 0;
  double Etmp2 = 0;

  for (int l = 0; l < N; l++)
  {
    double tmp_x = 0;
    double tmp_y = 0;
    double tmp_z = 0;

    for (int i = 0; i < N; i++)
    {
      if(l != i)
      {
        point_allias = pos[i];
        compute_allias(pos[l], point_allias);
        double o_ril = 1./r_ij(pos, point_allias, l);
        Etmp += gsl_pow_7(o_ril);

        tmp_x += (point_allias.x - pos[l].x) * gsl_pow_7(o_ril);
        tmp_y += (point_allias.y - pos[l].y) * gsl_pow_7(o_ril);
        tmp_z += (point_allias.z - pos[l].z) * gsl_pow_7(o_ril);
      }
    }
    Etmp2 += tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z;
  }
  return h2_2m * (-25/4 * b10 * Etmp2 + 10 * b5 * Etmp);
}

double E_kin2()
{
  double Etmp = 0;

  for (int l = 0; l < N; l++)
  {
    for (int i = 0; i < N; i++)
    if(i != l)
    {
      point_allias = pos[i];
      compute_allias(pos[l], point_allias);
      Etmp += gsl_pow_7(1/r_ij(pos, point_allias, l));
    }
  }
  return h2_2m * 5 * b5 * Etmp;
}

double E_pot()
{
  double Etmp = 0;

  for (int i = 0; i < N; i++)
  for (int j = 0; j < i; j++)
  {
    point_allias = pos[i];
    compute_allias(pos[j], point_allias);
    double rij6 = gsl_pow_6(1/r_ij(pos, point_allias, j));
    Etmp += 4*(rij6*rij6 - rij6);
  }

  return Etmp;
}

double E_pot2()
{
  double Etmp = 0;

  for (int i = 0;   i < N; i++)
  for (int j = i+1; j < N; j++)
  {
    point_allias = pos[i];
    compute_allias(pos[j], point_allias);
    double rij6 = gsl_pow_6(1/r_ij(pos, point_allias, j));
    Etmp += 4*(rij6*rij6 - rij6);
  }

  return Etmp;
}

double eloc;
double eloc2;
double var;
double var2;

double bMax = 1.3; //2.5;//1.25
double bMin = 1.1; //0.5;//1.05
double bStep = 0.01;

const double rho_exp = 21.86;
double rho_Min = rho_exp * 0.9;
double rho_Max = rho_exp * 1.1;
double n_rho = 10;

char filename[20] = "dati/initialPos.csv";
char filename1[20] = "dati/newPos.csv";
char filename2[40];
char filename3[40]; // = "dati/result.csv";
//char filename4[40];

int main()
{
  remove("dati/*.csv");
  remove("dati/particlePos/*.csv");

  for(double rho = rho_Min; rho <= rho_Max; rho += (rho_Max - rho_Min)/n_rho)
  {
    a = cbrt(32/rho) / 0.2556 / 2; // 4.442278245334543/2; experimental value
    printf("\n\n\na = %lg \nrho = %lg\n\n\n",a, rho);

    Delta = a * 0.05; //0.25

    placeParticles();

    FILE *fp;
    FILE *fp2;

    sprintf(filename3, "dati/results_rho%.5lg.csv", rho);
    fp2=fopen(filename3,"w+");

    int cnt = 0;

    for (b = bMin; b < bMax; b+= bStep)
    {
      cnt++;

      clock_t begin = clock();

      b5 = pow(b, 5);
      b10 = b5*b5;

      Eloc = 0;
      Eloc2 = 0;
      Var = 0;
      Var2 = 0;

      Nsamples = 0;
      tot = 0;

      sprintf(filename2, "dati/E%.4lg_rho%.5lg.csv", b, rho);
      fp=fopen(filename2,"w+");

      for (int t = 0; t < 100000; t++)
      {

        MC_Move();

        //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
        if(t > 2000)
        {
          tot += MC_acc();
          Nsamples++;
          pot = E_pot();
          pot2 = E_pot2();
          kin = E_kin();
          kin2 = E_kin2();
          Eloc += kin + pot; //calculate the energy
          Var += gsl_pow_2(kin + pot); //calculate the energy
          Eloc2 += kin2 + pot;
          Var2 += gsl_pow_2(kin2 + pot);

          /*
          if ((pot + kin)/Nsamples > 0.2)
          {
            printf("t: %d     ", t);
            printf("pot: %lg      kin: %lg\n", pot,kin);
            sprintf(filename4, "dati/particlePos/d%.4lg_t%d.csv", b,t);
            writeCSVparticles(filename4,pos);
          }

          if ((pot + kin2)/Nsamples > 0.2)
          {
            printf("t: %d     ", t);
            printf("pot: %lg      kin2: %lg\n", pot,kin2);
            sprintf(filename4, "dati/particlePos/d%.4lg_t%d.csv", b,t);
            writeCSVparticles(filename4,pos);
          }
          */

          eloc = Eloc/Nsamples;
          eloc2 = Eloc2/Nsamples;
          var = Var/Nsamples;
          var2 = Var2/Nsamples;

          //print data into the file
          fprintf(fp,"%ld %lg %lg %lg %lg %lg %lg %lg %lg \n",Nsamples, eloc, eloc2, sqrt((var - eloc*eloc)), sqrt((var2 - eloc2*eloc2)), kin, kin2, pot, pot2);
        }
        else
        {
          MC_acc();
        }



        if(t % 1000 == 999) printf("%d\n", t/1000);
      }

      clock_t end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

      printf("\n %lg \t %lg \n", eloc, sqrt((var - eloc*eloc)));
      printf("%lg \t %lg \n", eloc2, sqrt((var2 - eloc2*eloc2)));
      printf("%lg %lg\n", (double) tot/(double) Nsamples, b);
      printf("time = %lg s\n", time_spent);
      printf("%d/%d completati\n", cnt, (int) ((bMax - bMin)/bStep)+1);

      fprintf(fp2, "%lg %lg %lg %lg %lg %lg \n", b, eloc, eloc2, sqrt((var - eloc*eloc)), sqrt((var2 - eloc2*eloc2)), tot/(double) Nsamples);

      fclose(fp);
    }

    fclose(fp2);
  }
  return 0;
}
