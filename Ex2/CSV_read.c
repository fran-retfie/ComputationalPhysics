#include <stdarg.h> //va_*
#include <stdlib.h> //calloc
#include <string.h> // strlen, strcpy
#include <stdio.h>

// char* concat(int count, ...)
// {
//   va_list ap;
//   int i;
//
//   //Find required length to store merged string
//   int len = 1;
//   va_start(ap,count);
//   for(i=0; i<count; i++)
//   {
//     len += strlen(va_arg(ap,char*));
//   }
//   va_end(ap);
//
//   //Allocate memory to concat strings
//   char *merged = calloc(sizeof(char),len);
//   int null_pos = 0;
//
//   //Actually concatenate strings
//   va_start(ap,count);
//   for(i=0;i<count;i++)
//   {
//     char *s = va_arg(ap, char*);
//     strcpy(merged+null_pos, s);
//     null_pos += strlen(s);
//   }
//   va_end(ap);
//
//   return merged;
// }
//
// int readCSV(int row, int clm, float out[row][clm], char *filename)
// {
//   // build string based on columns of file
//   char *tmp;
//   tmp = concat(1, " %d,");
//   for(int i=0;i<clm-1;i++)
//   {
//     tmp = concat(2, tmp, "%f ,");
//   }
//   tmp = concat(2, tmp, "%f ");
//   printf("%s\n", tmp);
//
//   // Read .csv and fill the nxm matrix
//   FILE *fp;
//   fp = fopen(filename,"r");
//   int f,k;
//   for(int i=0;i<row;i++)
//   {
//     f = fscanf(fp, tmp, &k,&out[i][0],&out[i][1],&out[i][2],&out[i][3]);
//   }
//   return 0;
// }

int main()
{
  int n=4;
  int m=4;
  float Ex[n][m];


  char filename[12] = "Ex_terms.csv";




  readCSV(n,m,Ex,filename);
  for(int i=0;i<n;i++)
  {
    printf("%f ,%f ,%f ,%f\n", Ex[i][0],Ex[i][1],Ex[i][2],Ex[i][3]);
  }


  return 0;
}
