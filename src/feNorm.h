#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

class feNorm{

protected:
  feSpace *_intSpace;
  int cncGeoTag;
  feSpace *geoSpace;
  int nElmGeo;
  int dim;
  int nNodePerElem;
  std::vector<double> geoCoord;

  int _nQuad;
  int deg;
  std::vector<double> w;
  // std::vector<double> x;
  // std::vector<double> y;
  // std::vector<double> z;

  feFunction *_solRef;

  double norm;
public:
	feNorm(feSpace *intSpace, feMesh *mesh, int nQuad, feFunction *solRef = nullptr);
	~feNorm() {}

  double getNorm(){ return norm; }
  void computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
};

#endif