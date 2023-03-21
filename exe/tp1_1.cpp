
#include "feAPI.h"

//
// Analytical solution for the 2D diffusion problem
// defined on [0,1] x [0,1] : 
//
//  nabla^2 u(x,y) = 0
//          u(x,y) = g on the boundary
//
// where g = 0 on the left, right and bottom side
//   and g = x(1-x) on the top side 
//
double analyticalSolution(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double x = pos[0];
  double y = pos[1];

  double L = 1.;
  double H = 1.;

  double pi = M_PI;

  int nModes = 200;
  double res = 0.;
  double An, kn;
  for(int n = 1; n < nModes; n += 2) {
    kn = n * pi/L;
    An = 8. / (sinh(n*pi*H/L) * n*n*n * pi*pi*pi);
    res += An * sin(kn * x) * sinh(kn * y);
  }
  return res;
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);
  setVerbose(1);

  // ...

  petscFinalize();
  return 0;
}