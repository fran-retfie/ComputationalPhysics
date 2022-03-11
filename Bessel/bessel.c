#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define PI 3.141592654

double x,j,jj,n,nn,xmax,dx;
long int l;
FILE *bessj,*bessn;
char bj[20],bn[20];

void besselj()
{
  long il;
  double jp, jm, jz, norm;

  if(l==0 && x==0)
  {
    j=1;
    jj=0;
    return;
  }

  if(l>0 && x==0)
  {
    j=0;
    jj=0;
    return;
  }

  jm=0.;
  jz=9.9999e-300;
  il=230;

  while (il-- > 1)
  {
    jp=(((2.*il+1)*jz)/x)-jm;
    if((il-1)==l) j=jp;
    if((il-1)==l+1) jj=jp;
    jm=jz;
    jz=jp;
  }



  norm=sin(x)/x/jz;
  j=norm*j;
  jj=norm*jj;
  return;
}

void besseln()
{
  long il;
  double np, nm, nz;
  nm=(-cos(x))/x;
  nz=((-cos(x))/(x*x))-(sin(x))/x;

  if(l==0)
  {
    n=nm;
    nn=nz;
    return;
  }

  if(l==1)
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

void besseljteo()
{
  if(l==0) j = sin(x)/x;
  if(l==1) j = sin(x)/(x*x) - cos(x)/x;
  if(l==2) j = (3/(x*x)-1)*sin(x)/x - 3*cos(x)/(x*x);
  if(l==3) j = (15/(x*x*x)-6/x)*sin(x)/x - (15/(x*x)-1)*cos(x)/x;
  if(l==4) j = (1/x - 45/(x*x*x) +105/(x*x*x*x*x))*sin(x) + 5*(2/(x*x) - 21/(x*x*x*x))*cos(x);
  return;
}

void besselnteo()
{
  if(l==0) n = -cos(x)/x;
  if(l==1) n = -cos(x)/(x*x) - sin(x)/x;
  if(l==2) n = -(3/(x*x)-1)*cos(x)/x - 3*sin(x)/(x*x);
  if(l==3) n = (-15/(x*x*x)+6/x)*cos(x)/x - (15/(x*x)-1)*sin(x)/x;
  if(l==4) n = (-1/x + 45/(x*x*x) -105/(x*x*x*x*x))*cos(x) + 5*(2/(x*x) - 21/(x*x*x*x))*sin(x);
  return;
}

int main()
{
  printf("j function file name = ");
  scanf("%s",bj);
  printf("n function file name = ");
  scanf("%s",bn);

  bessj=fopen(bj,"w");
  bessn=fopen(bn,"w");
  x=0.1;
  xmax=15.0;
  dx=0.01;

  while(x<=xmax)
  {
    fprintf(bessj,"%lg\t",x);
    fprintf(bessn,"%lg\t",x);

    for(l=0;l<5;l++)
    {
      besseljteo();
      besselnteo();
      fprintf(bessj,"%lg\t",j);
      fprintf(bessn,"%lg\t",n);
    }
    x=x+dx;

    fprintf(bessj,"\n");
    fprintf(bessn,"\n");
  }
  fclose(bessn);
  fclose(bessj);
  return(0);
}
