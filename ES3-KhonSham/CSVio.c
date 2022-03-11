
//-------CSVio---------------------------
// A library for input output CVS files
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#include "CSVio.h"
#include <stdio.h>
#include <quadmath.h>

int readCSV(void)
{
  return 0;
}

int writeCSVint(char *filename, int *vec, int len1, int len2)
{
  FILE *fp;
  fp=fopen(filename,"w+");

  for(int i=0;i<len2;i++)
  {
      fprintf(fp,"\n%d",i+1);
      for(int j=0;j<len1;j++) fprintf(fp,",%i ", *(vec + len2*j + i) );
  }

  fclose(fp);
  return 0;
}

int writeCSVdouble(char *filename, double *vec, int len1, int len2)
{
  FILE *fp;
  fp=fopen(filename,"w+");

  for(int i=0;i<len2;i++)
  {
      fprintf(fp,"\n%d",i+1);
      for(int j=0;j<len1;j++) fprintf(fp,",%lg ", *(vec + len2*j + i) );
  }

  fclose(fp);
  return 0;
}

int writeCSVdouble_t(char *filename, double_t *vec, int len1, int len2, char *title)
{
  FILE *fp;
  fp=fopen(filename,"w+");

  fprintf(fp,"%s",title);

  for(int i=0;i<len2;i++)
  {
      fprintf(fp,"\n%d",i+1);
      for(int j=0;j<len1;j++)
      {
          if(sizeof(double_t) == sizeof(double))
            fprintf(fp,",%lg ", *(vec + len2*j + i) );
          else
            fprintf(fp,",%Qe ", *(vec + len2*j + i) );
      }
  }

  fclose(fp);
  return 0;
}
