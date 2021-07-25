
//-------DiffInt---------------------------
// A library for integrating differential equations
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#ifndef DiffInt
#define DiffInt

#include "doubleDef.h"

double RK4int(double t, double x, double y, double h, double (*f)(double, double, double));

double NumerovInt(double t, double y0, double y1, double h, double Ex, double (*f)(double,double));

double_t NumerovInt_t(int k, double_t y0, double_t y1, double_t h, double_t Ex, double_t (*f)(int,double_t));

double EulerInt(double t, double x, double y, double h, double (*f)(double, double, double));

#endif
