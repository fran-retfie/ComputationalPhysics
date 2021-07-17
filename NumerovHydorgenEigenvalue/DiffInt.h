
//-------DiffInt---------------------------
// A library for integrating differential equations
// Francesco Osti
// 06/03/2020
//-----------------------------------------

#ifndef DiffInt
#define DiffInt

#ifndef long_double_t
  typedef long double long_double_t; //__float128  long_double_t;
#endif

double RK4int(double t, double x, double y, double h, double (*f)(double, double, double));

double NumerovInt(double t, double y0, double y1, double h, double Ex, double (*f)(double,double));

long_double_t NumerovIntLD(long_double_t t, long_double_t y0, long_double_t y1, long_double_t h, long_double_t Ex, long_double_t (*f)(long_double_t,long_double_t));

double EulerInt(double t, double x, double y, double h, double (*f)(double, double, double));

#endif
