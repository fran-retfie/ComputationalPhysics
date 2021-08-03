
//-------CSVio---------------------------
// A library for input output CVS files
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#include <stdarg.h> //va_*
#include <stdlib.h> //calloc
#include <string.h> // strlen, strcpy
#include <stdio.h>

char* concat(int count, ...)
{
  va_list ap;
  int i;

  //Find required length to store merged string
  int len = 1;
  va_start(ap,count);
  for(i=0; i<count; i++)
  {
    len += strlen(va_arg(ap,char*));
  }
  va_end(ap);

  //Allocate memory to concat strings
  char *merged = calloc(sizeof(char),len);
  int null_pos = 0;

  //Actually concatenate strings
  va_start(ap,count);
  for(i=0;i<count;i++)
  {
    char *s = va_arg(ap, char*);
    strcpy(merged+null_pos, s);
    null_pos += strlen(s);
  }
  va_end(ap);

  return merged;
}

int readCSV(int row, int clm, float out[row][clm], char *filename)
{
  // build string based on columns of file
  char *tmp;
  tmp = concat(1, " %d,");
  for(int i=0;i<clm-1;i++)
  {
    tmp = concat(2, tmp, "%f ,");
  }
  tmp = concat(2, tmp, "%f ");
  // printf("%s\n", tmp);

  // Read .csv and fill the nxm matrix
  FILE *fp;
  fp = fopen(filename,"r");
  int f,k;
  for(int i=0;i<row;i++)
  {
    f = fscanf(fp, tmp, &k,&out[i][0],&out[i][1],&out[i][2],&out[i][3]);
  }
  fclose(fp);
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
