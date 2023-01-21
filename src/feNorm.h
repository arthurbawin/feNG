#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

typedef enum
{
  L1,
  L1_ERROR,
  L2,
  L2_ERROR,
  LINF,
  LINF_ERROR,
  SEMI_H1,
  SEMI_H1_ERROR,
  H1,
  H1_ERROR,
  INTEGRAL,
  DOT_PRODUCT,
  AREA
} normType;

class feNorm;

// Create a norm object and perform safety checks.
// This is the recommended way of creating an FE space.
// Call within feCheck( ... ) to safely exit if a problem arises.
feStatus createNorm(feNorm *&norm,
  normType type,
  const std::vector<feSpace*> &spaces,
  feSolution *sol,
  feFunction *scalarSolution,
  feVectorFunction *vectorSolution = nullptr);

//
// Computes a variety of integrals over a geometric connectivity
//
class feNorm
{
protected:
  normType _type;

  std::vector<feSpace *> _spaces;

  const feSolution *_solution;

  const feCncGeo *_cnc;
  int _nElm;

  std::vector<double> _w;
  int _nQuad;
  const std::vector<double> &_J;

  std::vector<std::vector<feInt>> _adr;
  std::vector<std::vector<double>> _localSol;

  std::vector<double> _pos;
  std::vector<double> _geoCoord;

  feSpace *_geoSpace;
  feFunction *_scalarSolution;
  feVectorFunction *_vectorSolution = nullptr;

public:
  feNorm(normType type, const std::vector<feSpace*> &spaces, feSolution *sol, feFunction *solRef = nullptr,
         feVectorFunction *VecSolRef = nullptr);
  ~feNorm() {}

  double compute();
  double compute(normType type);

private:
  void initializeLocalSolutionOnSpace(int iSpace, int iElm);

  double computeLpNorm(int p, bool error = false);
  double computeL1Norm(bool error = false);
  double computeL2Norm(bool error = false);
  double computeLInfNorm(bool error = false);
  double computeH1Norm(bool error = false);

  double computeArea();
  double computeIntegral();
  double computeIntegralDotProduct();

  feStatus computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                feMesh *mesh, feMetaNumber *refMN,
                                                feSolution *refSol, feMesh *refMesh,
                                                const std::vector<feSpace *> refSpaces,
                                                double &res);
};

#endif
