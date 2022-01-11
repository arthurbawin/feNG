#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

class feNorm {
protected:
  feSpace *_intSpace;

  int cncGeoTag;
  feSpace *geoSpace;
  int nElmGeo;
  int dim;
  int nNodePerElem;
  std::vector<double> geoCoord;

  int _nQuad;
  int _degQuad;
  int deg;
  std::vector<double> w;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  feFunction *_solRef;
  feVectorFunction *_solRefVec;

  double norm;

public:
  feNorm(feSpace *intSpace, feMesh *mesh, int degQuad, feFunction *solRef = nullptr);

  ~feNorm() {}

  double getNorm() { return norm; }
  void computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeL2Norm0D(feSolution *sol);
  void computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  void computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, 
    feMetaNumber *refMN, feSolution *refSol, feMesh *refMesh, const std::vector<feSpace *> refSpaces);
};

#endif
