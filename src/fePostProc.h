#ifndef _FEPOSTPROC_
#define _FEPOSTPROC_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

class fePostProc
{
protected:
  feSpace *_intSpace;
  feSpace *_geoSpace;

  feFunction *_referenceSolution;
  feVectorFunction *_referenceVectorSolution;

  feMesh *_mesh;
  feMetaNumber *_metaNumber;
  feCncGeo *_cnc;
  std::vector<int> _adr;
  std::vector<double> _sol;

public:
  fePostProc(feSpace *intSpace, feMesh *mesh, feMetaNumber *metaNumber,
             feFunction *referenceSolution = nullptr,
             feVectorFunction *referenceVectorSolution = nullptr)
    : _intSpace(intSpace), _geoSpace(intSpace->getCncGeo()->getFeSpace()),
      _referenceSolution(referenceSolution), _referenceVectorSolution(referenceVectorSolution),
      _mesh(mesh), _metaNumber(metaNumber), _cnc(intSpace->getCncGeo())
  {
    _adr.resize(_intSpace->getNbFunctions());
    _sol.resize(_adr.size());
  };

  ~fePostProc() {}

  double computeMeasure();
  double computeSolutionIntegral(feSolution *sol);
  double computeFunctionIntegral(feFunction *fun, double t);

  double computeL2ErrorNorm(feSolution *sol);

  // void computeL2Norm0D(feSolution *sol);
  // void computeIntegral(feMetaNumber *metaNumber, feSolution *sol, feMesh *mesh, feFunction *fun);
  // void computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol, feMesh
  // *mesh,
  //   feMetaNumber *refMN, feSolution *refSol, feMesh *refMesh, const std::vector<feSpace *>
  //   refSpaces);
};

#endif
