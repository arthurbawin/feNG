#include "feTriangle.h"

void Triangle::xyz2uvw(double xyz[3], double uvw[3])
{
  const double O[3] = {_v[0]->x(), _v[0]->y(), _v[0]->z()};

  const double d[3] = {xyz[0] - O[0], xyz[1] - O[1], xyz[2] - O[2]};
  const double d1[3] = {_v[1]->x() - O[0], _v[1]->y() - O[1], _v[1]->z() - O[2]};
  const double d2[3] = {_v[2]->x() - O[0], _v[2]->y() - O[1], _v[2]->z() - O[2]};

  const double Jxy = d1[0] * d2[1] - d1[1] * d2[0];
  const double Jxz = d1[0] * d2[2] - d1[2] * d2[0];
  const double Jyz = d1[1] * d2[2] - d1[2] * d2[1];

  if((std::abs(Jxy) > std::abs(Jxz)) && (std::abs(Jxy) > std::abs(Jyz))) {
    uvw[0] = (d[0] * d2[1] - d[1] * d2[0]) / Jxy;
    uvw[1] = (d[1] * d1[0] - d[0] * d1[1]) / Jxy;
  } else if(std::abs(Jxz) > std::abs(Jyz)) {
    uvw[0] = (d[0] * d2[2] - d[2] * d2[0]) / Jxz;
    uvw[1] = (d[2] * d1[0] - d[0] * d1[2]) / Jxz;
  } else {
    uvw[0] = (d[1] * d2[2] - d[2] * d2[1]) / Jyz;
    uvw[1] = (d[2] * d1[1] - d[1] * d1[2]) / Jyz;
  }
  uvw[2] = 0.0;
}