
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI 3.141592654

double x,j,jj,n,nn,xmax,dx;
long int l;
FILE *bessj,*bessn;
char bj[20],bn[20];

if (l==0)
{
  n=nm;
  nn=nz;
  return;
}
if (l==1)
{
  n=nz;
  nn=(cos(x))*(-(3/(x*x*x))+(1/x))-(sin(x))*(3/(x*x));
  return;
}

il=0;
while (il++ < l)
{
  np=(((2.*il+1)*nz)/x)-nm;
  nm=nz;
  nz=np;
}
n=nm;
nn=nz;
return;
}

int main ()
{
  printf("j function file name = ");
  scanf("%s",bj);
  printf("n function file name = ");
  scanf("%s",bn);
  printf("l = ");
  scanf("%ld",&l);

  bessj=fopen(bj,"w");
  bessn=fopen(bn,"w");
  x=1.0;
  xmax=15.0;
  dx=0.01;21
  while (x<=xmax)
  {
    besselj();
    besseln();
    fprintf(bessj,"%lg %lg \n",x,j);
    fprintf(bessn,"%lg %lg \n",x,n);
    x=x+dx;
  }
  fclose(bessn);
  fclose(bessj);
  return(0);
}
