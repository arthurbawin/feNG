
#include "feAPI.h"
#include "../tests.h"

#include <gtest/gtest.h>

double fxSol(const feFunctionArguments &args, const std::vector<double> &par)
{
  const double x = args.pos[0];
  const int deg = par[0];
  return pow((x+1), deg);
  // return pow(x, deg);
}

double fySol(const feFunctionArguments &args, const std::vector<double> &par)
{
  const double y = args.pos[1];
  const int deg = par[0];
  return pow((y+1), deg);
  // return pow(y, deg);
}

int computeQuadrature()
{
  for(int ideg = 0; ideg <= 25; ++ideg)
  {
    // Precision decreases for high-order polynomials
    const double err_tolerance = (ideg < 12) ? 1e-13 : 5e-10;

    feFunction funx(fxSol, {(double) ideg});
    feFunction funy(fySol, {(double) ideg});

    // Define a quadrature rule to integrate polynomials of order up to 2 on reference triangle
    feQuadrature rule(ideg, geometryType::TRI);
    std::vector<double> wQ = rule.getWeights();
    std::vector<double> xsiQ = rule.getXPoints();
    std::vector<double> etaQ = rule.getYPoints();

    feFunctionArguments args;

    double Ix = 0., Iy = 0.;
    for(size_t k = 0; k < wQ.size(); ++k) {
      args.pos[0] = xsiQ[k];
      args.pos[1] = etaQ[k];
      Ix += funx(args) * wQ[k];
      Iy += funy(args) * wQ[k];
    }

    const double p = (double) ideg;
    double Iexact = (pow(2., ideg+2) - p - 3.) / (p*p + 3.*p + 2.);
    // double Iexact = 1./((p+1.)*(p+2.));
    double errx = fabs(Ix - Iexact);
    double erry = fabs(Iy - Iexact);

    // printf("deg = %2d : error on funx = %+-1.15e - error on funy = %+-1.15e\n", ideg, errx, erry);
    if(errx > err_tolerance) return 1;
    if(erry > err_tolerance) return 1;
  }
  return 0;
}

TEST(Quadrature, QuadratureTriangle) {
  ASSERT_TRUE(computeQuadrature() == 0);  
}

