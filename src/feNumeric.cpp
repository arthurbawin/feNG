//
// Adapted from Gmsh's Numeric.cpp
//

#include "feMessage.h"
#include "feNumeric.h"

double det2x2(double mat[2][2])
{
  return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
}

double det3x3(double mat[3][3])
{
  return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
          mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
          mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));
}

double inv2x2(double mat[2][2], double inv[2][2])
{
  const double det = det2x2(mat);
  if(det) {
    double ud = 1. / det;
    inv[0][0] = mat[1][1] * ud;
    inv[1][0] = -mat[1][0] * ud;
    inv[0][1] = -mat[0][1] * ud;
    inv[1][1] = mat[0][0] * ud;
  }
  else {
    feErrorMsg(FE_STATUS_ERROR, "Singular matrix 2x2");
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++) inv[i][j] = 0.;
  }
  return det;
}

double inv3x3(double mat[3][3], double inv[3][3])
{
  double det = det3x3(mat);
  if(det) {
    double ud = 1. / det;
    inv[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) * ud;
    inv[1][0] = -(mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) * ud;
    inv[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) * ud;
    inv[0][1] = -(mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) * ud;
    inv[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * ud;
    inv[2][1] = -(mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) * ud;
    inv[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * ud;
    inv[1][2] = -(mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]) * ud;
    inv[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) * ud;
  }
  else {
    feErrorMsg(FE_STATUS_ERROR, "Singular matrix 3x3");
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++) inv[i][j] = 0.;
  }
  return det;
}