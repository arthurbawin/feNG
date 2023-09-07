#ifndef _FENORMV2_
#define _FENORMV2_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feQuadrature.h"
#include "feSolutionContainerV2.h"

class feSolutionContainerV2;

class feNormV2
{
protected:
  std::vector<feSpace *> _intSpaces;
  feSpace *_geoSpace;
  feMesh *_mesh;
  feMetaNumber *_metaNumber;
  feCncGeo *_cnc;

  std::string _normType;
  std::vector<feFunction *> _referenceSolutions;

  int _nElm;
  int _nNodePerElem;

  int _nQuad;
  std::vector<double> &_wQuad;
   

public:
  feNormV2(std::vector<feSpace *> spaces, feMesh *mesh, feMetaNumber *metaNumber, std::string normType, std::vector<feFunction *> referenceSolutions = {nullptr}):
      _intSpaces(spaces), _mesh(mesh), _metaNumber(metaNumber), _normType(normType), _referenceSolutions(referenceSolutions), _nElm(_intSpaces[0]->getNbElm())
      , _nNodePerElem(_intSpaces[0]->getNbNodePerElem()), _geoSpace(_intSpaces[0]->getCncGeo()->getFeSpace()), _cnc(_intSpaces[0]->getCncGeo())
      , _nQuad(_geoSpace->getNbQuadPoints()), _wQuad(_intSpaces[0]->getQuadratureWeights()) {};

  double compute(double tC, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution=nullptr);

  double computeL2ErrorNorm(double tC, feSolutionContainerV2 *solutionContainer);
  double computeL2ErrorNormVec(double t, feSolutionContainerV2 *solutionContainer);
  double computeL2ErrorNormFromMaxDC(double t, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution);
  double computeL2ErrorNormVecFromMaxDC(double t, feSolutionContainerV2 *solutionContainer, feSolutionContainerV2 *referenceSolution);

  double integrationLagrangeMupltiplicator(double t, feSolutionContainerV2 *solutionContainer);



  ~feNormV2() {}
};




#endif