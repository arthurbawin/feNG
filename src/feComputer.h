#ifndef _FECOMPUTER_
#define _FECOMPUTER_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

class feComputer
{
protected:
  feSpace *_intSpace;
  int _nNodePerElem;
  std::vector<feSpace *> _VecfeSpace;
  feSpace *_geoSpace;

  feFunction *_referenceSolution;
  feVectorFunction *_referenceVectorSolution;

  feMesh *_mesh;
  feMetaNumber *_metaNumber;
  feCncGeo *_cnc;

  std::string _computerType;

public:
  std::vector<double> _Result;

  feComputer(feSpace *intSpace, feMesh *mesh, feMetaNumber *metaNumber,
             std::string computerType = "", feFunction *referenceSolution = nullptr,
             feVectorFunction *referenceVectorSolution = nullptr)
    : _intSpace(intSpace), _geoSpace(intSpace->getCncGeo()->getFeSpace()),
      _referenceSolution(referenceSolution), _computerType(computerType),
      _referenceVectorSolution(referenceVectorSolution), _mesh(mesh), _metaNumber(metaNumber),
      _cnc(intSpace->getCncGeo()), _nNodePerElem(_intSpace->getNbNodePerElem()){};

  feComputer(std::vector<feSpace *> &VecfeSpace, feMesh *mesh, feMetaNumber *metaNumber,
             std::string computerType = "", feFunction *referenceSolution = nullptr,
             feVectorFunction *referenceVectorSolution = nullptr)
    : _VecfeSpace(VecfeSpace), _geoSpace(VecfeSpace[0]->getCncGeo()->getFeSpace()),
      _referenceSolution(referenceSolution), _computerType(computerType),
      _referenceVectorSolution(referenceVectorSolution), _mesh(mesh), _metaNumber(metaNumber),
      _cnc(VecfeSpace[0]->getCncGeo()), _nNodePerElem(VecfeSpace[0]->getNbNodePerElem()){};

  ~feComputer() {}

  void Compute(feSolution *sol, int i);

  std::vector<double> &getResult() { return _Result; };

  // feCncGeo *getCncGeo(){return _cnc; };
  feSpace *getIntSpace(){return _intSpace;};
  double computeL2ErrorNorm(feSolution *sol);
  double computeL2ErrorNormVec(feSolution *sol);
  double computeH1ErrorNormVec(feSolution *sol);
  double computeSemiH1ErrorNormVec(feSolution *sol);
  double computeNormL2Lambda(feSolution *sol);
  double computeIntLagrangeMultiplier(feSolution *sol);
  double computeIntMultiplierNS(feSolution *sol);
};

#endif