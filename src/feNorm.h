#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

class feNorm
{
protected:
  std::vector<feSpace *> _VecfeSpace;
  feSpace *_intSpace; // assossiated with _VecfeSpace[0]
  int _NbFields = 1;
  int _cncGeoTag;
  feSpace *_geoSpace;
  int _nElmGeo;
  int _dim;
  int _nNodePerElem;
  std::vector<double> _geoCoord;

  int _nQuad;
  int _degQuad;
  int _deg;
  std::vector<double> _w;
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _z;

  feFunction *_solRef;
  feVectorFunction *_VecSolRef = nullptr;
  double _norm;
  std::string _TypeNorm;

public:
  feNorm(std::vector<feSpace *> &VecfeSpace, feMesh *mesh, int degQuad,
         feFunction *solRef = nullptr, feVectorFunction *VecSolRef = nullptr,
         std::string TypeNorm = "");
  feNorm(feSpace *intSpace, feMesh *mesh, int degQuad, feFunction *solRef = nullptr,
         feVectorFunction *VecSolRef = nullptr, std::string TypeNorm = "");
  ~feNorm() {}
  // ~feNormVec() {}

  double getNorm() { return _norm; }
  int getNbFields() { return _NbFields; }
  std::string getTypeNorm() { return _TypeNorm; }
  void computeL2Norm(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeL2Norm_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, double scaling = 1.0);
  void computeL2NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeH1NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeSemiH1NormVec(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeNormLambda(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeNormLambdaNS(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntFluxVeloc(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntFluxNx(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntFluxNy(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntegralNum(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  // void computeIntFlux(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);

  void computeL2Norm0D(feSolution *sol);
  void computeArea(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh);
  void computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  void computeScalarProduct(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  void computeL2Norm_f_minus_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun, double scaling = 1.);

  feStatus computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                feMesh *mesh, feMetaNumber *refMN,
                                                feSolution *refSol, feMesh *refMesh,
                                                const std::vector<feSpace *> refSpaces);
};

#endif
