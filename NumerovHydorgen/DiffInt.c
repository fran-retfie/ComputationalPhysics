
//-------DiffInt---------------------------
// A library for integrating differential equations
// Francesco Osti
// 06/03/2020
//-----------------------------------------
#include "stdio.h"
#include "DiffInt.h"

double RK4int(double t, double x, double y, double h, double (*f)(double, double, double))
{
  double h2 = h*0.5;
  double k1 = (*f)(t, x, y);
  double k2 = (*f)(t + h2, x + k1*h2, y + k1*h2);
  double k3 = (*f)(t + h2, x + k2*h2, y + k2*h2);
  double k4 = (*f)(t + h, x + k3*h , y + k3*h);

  return y + h/6.0*(k1 + 2*k2 + 2*k3 + k4);
}

double NumerovInt(double x, double y0, double y1, double h, double Ex, double (*f)(double,double))
{
  double const o12 = 1/12.0;
  double h2 = h*h;
  return (2*(1-5*h2*f(x,Ex)*o12)*y0 - (1+h2*f(x-h,Ex)*o12)*y1)/(1+h2*f(x+h,Ex)*o12);
}

long_double_t NumerovIntLD(long_double_t x, long_double_t y0, long_double_t y1, long_double_t h, long_double_t Ex, long_double_t (*f)(long_double_t,long_double_t))
{
  long_double_t const o12 = 1/12.0;
  long_double_t h2 = h*h;
  return (2*(1-5*h2*f(x,Ex)*o12)*y0 - (1+h2*f(x-h,Ex)*o12)*y1)/(1+h2*f(x+h,Ex)*o12);
}

double EulerInt(double t, double x, double y, double h, double (*f)(double, double, double))
{
  return y + (*f)(t, x, y) * h;
}
