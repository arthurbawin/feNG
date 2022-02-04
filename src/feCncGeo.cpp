#include "feCncGeo.h"
#include "feSpace.h"
#include "feMesh.h"

void feCncGeo::computeJacobians()
{
  int nQuad = _space->getNbQuadPoints();
  _J.resize(_nElm * nQuad);

  std::vector<double> geoCoord;

  if(_space->getDim() == 1) {
    std::vector<double> j(3, 0.0);
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      geoCoord = _mesh->getCoord(_tag, iElm);
      for(int k = 0; k < nQuad; ++k) {
        // _space->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, j);
        // _J[nQuad * iElm + k] = j[0];
        std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
        _space->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        _J[nQuad * iElm + k] = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
      }
    }
  } else if(_space->getDim() == 2) {
    std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
    std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      geoCoord = _mesh->getCoord(_tag, iElm);
      for(int k = 0; k < nQuad; ++k) {
        _space->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        _space->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
        if(_J[nQuad * iElm + k] < 0) {
          printf("In feCncGeo::computeJacobians : Error - Element jacobian = %+-12.12e\n",
                 _J[nQuad * iElm + k]);
          exit(-1);
        }
      }
    }
  } else if(_space->getDim() == 0) {
    _J[0] = 1.0;
  }

  else {
    printf("TODO : Implement jacobian for 3D elements.\n");
  }
}