#ifndef _ELLIPSETOOLBOX_
#define _ELLIPSETOOLBOX_

#include <vector>
#include <iostream>

/* This is essentially a C++ translation of these matlab conic functions :

  Quadratic curves and quadric surfaces in implicit form
  Copyright 2010 Levente Hunyadi

  Author : Arthur Bawin
*/

// Sign function
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

/* p is the coefficients of the implicit form a*x^2 + b*x*y + c*y^2 + d*x + f*y + g = 0
   a and b are the explicit coefficients : x^2/a^2 + y^2/b^2 = 1
   R is the rotation matrix
   t is the translation vector */
bool getExplicitEllipse(std::vector<double> &p, double *a, double *b, std::vector<double> &R,
                        std::vector<double> &t);

bool getEllipsePoints(double a, double b, double c, double xC, double yC, std::vector<double> &x,
                      std::vector<double> &y);

#endif