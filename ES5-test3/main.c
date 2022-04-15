
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

char filename[20] = "dati/initialPos.csv";
char filename1[20] = "dati/newPos.csv";
char filename2[20];
char filename3[20] = "dati/result.csv";
char filename4[40];

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

long Nsamples;
long tot;



double const h2_2m = 0.09075588736790198; // 662.1821253339011;

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

// from volpx
int sign(double x)
{
  return x > 0 ? 1 : -1;
}

double computeAlias(struct point pos_i, struct point pos_j)//compute the alias particle wrt j
{
  if (fabs(pos_j.x - pos_i.x) > a) pos_i.x += sign(pos_j.x - pos_i.x)*(2*a);
  if (fabs(pos_j.y - pos_i.y) > a) pos_i.y += sign(pos_j.y - pos_i.y)*(2*a);
  if (fabs(pos_j.z - pos_i.z) > a) pos_i.z += sign(pos_j.z - pos_i.z)*(2*a);

  double dx = pos_j.x - pos_i.x;
  double dy = pos_j.y - pos_i.x;
  double dz = pos_j.z - pos_i.z;

  double r = sqrt(dx*dx + dy*dy + dz*dz);
  //printf("%lg\n",r < a ? r : -r);
  return r < a ? r : -r;
}

void apply_periodic_bounds()
{
  for(int i=0;i<N;i++)
  {
    posNew[i].x -= (2*a)*floor(posNew[i].x/(2*a));
    posNew[i].y -= (2*a)*floor(posNew[i].y/(2*a));
    posNew[i].z -= (2*a)*floor(posNew[i].z/(2*a));
  }
}

void MC_Move()
{
  for(int i=0;i<N;i++)
  {
    posNew[i].x = pos[i].x + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].y = pos[i].y + gsl_ran_flat(r, -Delta, Delta);
    posNew[i].z = pos[i].z + gsl_ran_flat(r, -Delta, Delta);


/*
    while(posNew[i].x > a)  posNew[i].x -= 2*a;
    while(posNew[i].y > a)  posNew[i].y -= 2*a;
    while(posNew[i].z > a)  posNew[i].z -= 2*a;
    while(posNew[i].x < -a) posNew[i].x += 2*a;
    while(posNew[i].y < -a) posNew[i].y += 2*a;
    while(posNew[i].z < -a) posNew[i].z += 2*a;
*/
  }
  apply_periodic_bounds(posNew);
}

double r_ij(struct point points[N], int i, int j)
{
  double dx = fabs(points[j].x - points[i].x);
  double dy = fabs(points[j].y - points[i].y);
  double dz = fabs(points[j].z - points[i].z);

  while(dx > a) dx = 2*a - dx;
  while(dy > a) dy = 2*a - dy;
  while(dz > a) dz = 2*a - dz;

  return gsl_hypot3(dx, dy, dz);
}

double normScalarProd(struct point points[N], int l, int i, int j)
{
  double dxi = points[l].x - points[i].x;
  double dyi = points[l].y - points[i].y;
  double dzi = points[l].z - points[i].z;

  while(dxi > a) dxi = dxi - 2*a;
  while(dyi > a) dyi = dyi - 2*a;
  while(dzi > a) dzi = dzi - 2*a;

  while(dxi < -a) dxi = 2*a + dxi;
  while(dyi < -a) dyi = 2*a + dyi;
  while(dzi < -a) dzi = 2*a + dzi;

  double dxj = points[l].x - points[j].x;
  double dyj = points[l].y - points[j].y;
  double dzj = points[l].z - points[j].z;

  while(dxj > a) dxj = dxj - 2*a;
  while(dyj > a) dyj = dyj - 2*a;
  while(dzj > a) dzj = dzj - 2*a;

  while(dxj < -a) dxj = 2*a + dxj;
  while(dyj < -a) dyj = 2*a + dyj;
  while(dzj < -a) dzj = 2*a + dzj;

  double prod = dxi * dxj + dyi * dyj + dzi * dzj;

  return prod/(r_ij(posNew, i, l) * r_ij(posNew, j, l));
}

int hardSphere(double r)
{
  if (r < 0.01) // distance normalized on sigma
  {
    return 1;
  }
  else return 0;

}

int MC_acc()
{
  double acc;
  double Orij;
  double Nrij;
  double sumPsi = 0;
  int check = 0; // check is > 0 if two particles are closer than 1 sigma (i.e. 1)

  for (int j = 0; j < N; j++)
  for (int i = 0; i < j; i++)
  {
    Orij = computeAlias(pos[i],pos[j]);//r_ij(pos, i, j);
    Nrij = computeAlias(posNew[i],posNew[j]);//r_ij(posNew, i, j);
    if (Nrij > 0) sumPsi += gsl_pow_5(1/Orij) - gsl_pow_5(1/Nrij);

    //check += hardSphere(Nrij);
  }

  if (check == 0)
  {
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
  else return 0;
}
/*

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
        double o_ril = 1./r_ij(pos, i, l);
        Etmp += gsl_pow_7(o_ril);

        for (int j = 0; j < N; j++)
        {
          if(l != j)
            Etmp2 += normScalarProd(pos, l, i, j) * gsl_pow_6(o_ril) * gsl_pow_6(1/r_ij(pos, j, l));
        }
      }
    }
  }
  return h2_2m * (-25/4 * b10 * Etmp2 + 10 * b5 * Etmp);
}
*/

double E_kin()
{
  double Etmp = 0;
  double sum = 0;
  double tmp_x=0;
  double tmp_y=0;
  double tmp_z=0;


  for (int l = 0; l < N; l++)
  {
    for (int i = 0; i < N; i++)
    {
      if(l != i)
      {
        double o_ril = computeAlias(pos[i],pos[l]);//r_ij(pos, i, l);
        if (o_ril>0)
        {
          Etmp += gsl_pow_7(1./o_ril);
          tmp_x += (pos[l].x - pos[i].x)/o_ril;
          tmp_y += (pos[l].y - pos[i].y)/o_ril;
          tmp_z += (pos[l].z - pos[i].z)/o_ril;
        }
      }
    }
  }
  return h2_2m * (-25/4 * b10 * (tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z) + 10 * b5 * Etmp);
}

double E_kin2()
{
  double Etmp = 0;
  double o_ril = 0;

  for (int l = 0; l < N; l++)
  {
    for (int i = 0; i < N; i++)
    if(i != l)
    {
      o_ril = computeAlias(pos[i],pos[l]);
      if (o_ril>0) Etmp += gsl_pow_7(1./o_ril);//r_ij(pos, i, l));
    }
  }
  return h2_2m * 5 * b5 * Etmp;
}

double E_pot()
{
  double Etmp = 0;
  double rij = 0;

  for (int i = 0; i < N; i++)
  for (int j = i+1; j < N; j++)
  {
    rij = computeAlias(pos[i],pos[j]);
    if (rij>0)
    {
      double rij6 = gsl_pow_int(rij, -6);// gsl_pow_int(r_ij(pos, i, j), -6);
      Etmp += 4*(rij6*rij6 - rij6);
    }
  }

  return Etmp;
}

double eloc;
double eloc2;
double var;
double var2;

double bMax = 1.4;//1.25
double bMin = 1.;//1.05
double bStep = 0.01;

int main()
{
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  a = 4.442278245334543/2.;
  Delta = a * 0.03;

  placeParticles();
  //apply_periodic_bounds();

  remove("dati/*.csv");
  remove("dati/particlePos/*.csv");

  FILE *fp;
  FILE *fp2;

  printf("%lg\n", rint((-a-0.1)/(2*a)));
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

    sprintf(filename2, "dati/E%.4lg.csv", b);
    fp=fopen(filename2,"w+");

    for (int t = 0; t < 12000; t++)
    {

      MC_Move();

      //wait until thermalization (t = 150, Delta = 0.002), see drawAcc.p
      if(t > 2000)
      {
        tot += MC_acc();
        Nsamples++;
        double pot = E_pot();
        double kin = E_kin();
        double kin2 = E_kin2();
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
        fprintf(fp,"%ld %lg %lg %lg %lg \n",Nsamples, eloc, eloc2, sqrt((var - eloc*eloc)), sqrt((var2 - eloc2*eloc2)));
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

  return 0;
}
