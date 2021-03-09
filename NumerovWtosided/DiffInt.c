
//-------DiffInt---------------------------
// A library for integrating differential equations
// Francesco Osti
// 06/03/2020
//-----------------------------------------
#include "stdio.h"
#include "DiffInt.h"

double RK4int(double t, double x, double y, double h, double (*f)(double, double, double))
{
  double k1 = (*f)(t, x, y);
  double k2 = (*f)(t + h/2.0, x + k1*h/2.0, y + k1*h/2.0);
  double k3 = (*f)(t + h/2.0, x + k2*h/2.0, y + k2*h/2.0);
  double k4 = (*f)(t + h, x + k3*h , y + k3*h);

  return y + h/6.0*(k1 + 2*k2 + 2*k3 + k4);
}

double NumerovInt(double x, double y0, double y1, double h, double Ex, double (*f)(double,double))
{
  double h2 = h*h;
  return (2*(1-5*h2*f(x,Ex)/12)*y0 - (1+h2*f(x-h,Ex)/12)*y1)/(1+h2*f(x+h,Ex)/12);
}

double NumerovIntBackward(double x, double y0, double y1, double h, double Ex, double (*f)(double,double))
{
  double h2 = h*h;
  return (2*(1-5*h2*f(x,Ex)/12)*y0 - (1+h2*f(x+h,Ex)/12)*y1)/(1+h2*f(x-h,Ex)/12);
}

double EulerInt(double t, double x, double y, double h, double (*f)(double, double, double))
{
  return y + (*f)(t, x, y) * h;
}
