
//-------DiffInt---------------------------
// A library for integrating differential equations
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#ifndef DiffInt
#define DiffInt

double RK4int(double t, double x, double y, double h, double (*f)(double, double, double));

double NumerovInt(double t, double y0, double y1, double h, double Ex, double (*f)(double,double));

double EulerInt(double t, double x, double y, double h, double (*f)(double, double, double));

#endif
