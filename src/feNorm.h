#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

typedef enum { L1, L2, LINF, H1 } normType;

//
// Computes a norm ||.||_normType over a geometric connectivity
//
class feNorm
{
protected:
  std::vector<feSpace *> _spaces;

  const feCncGeo *_cnc;
  int _nElm;
  int _nNodePerElem;

  std::vector<double> _w;
  int _nQuad;
  const std::vector<double> &_J;

  std::vector<std::vector<feInt>> _adr;
  std::vector<std::vector<double>> _sol;

  std::vector<double> _pos;
  std::vector<double> _geoCoord;

  feSpace *_geoSpace;
  feFunction *_scalarSolution;
  feVectorFunction *_vectorSolution = nullptr;

  std::vector<double> _normHistory;

  double _norm;

public:
  feNorm(const std::vector<feSpace*> &spaces, feMesh *mesh, int degQuad, feFunction *solRef = nullptr,
         feVectorFunction *VecSolRef = nullptr);
  ~feNorm() {}

  double getNorm() { return _norm; }
  void computeL2Norm(feSolution *sol, feMesh *mesh);
  void computeL2Norm_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, double scaling = 1.0);

  void computeL2Norm0D(feSolution *sol);
  void computeArea();
  void computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  void computeScalarProduct(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  void computeL2Norm_f_minus_uh(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun, double scaling = 1.);

  feStatus computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                feMesh *mesh, feMetaNumber *refMN,
                                                feSolution *refSol, feMesh *refMesh,
                                                const std::vector<feSpace *> refSpaces);
};

#endif
