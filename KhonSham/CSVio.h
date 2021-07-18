

//-------CSVio---------------------------
// A library for input output CVS files
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#ifndef CVSio
#define CVSio

#include "doubleDef.h"

int readCSV(void);

int writeCSVint(char *filename ,int *vec, int len1, int len2);
int writeCSVdouble(char *filename ,double *vec, int len1, int len2);
int writeCSVdouble_t(char *filename ,double_t *vec, int len1, int len2);

#endif
