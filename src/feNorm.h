#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

typedef enum
{
  // In the following, cnc refers to the geometric connectivity
  // over which the norm or integral is computed.

  // Compute ||uh||_L1(cnc)
  L1,
  // Compute ||u - uh||_L1(cnc), with u a user-defined scalar field
  L1_ERROR,
  // Compute ||uh||_L2(cnc)
  L2,
  VECTOR_L2,
  // Compute ||u - uh||_L1(cnc), with u a user-defined scalar field
  L2_ERROR,
  VECTOR_L2_ERROR,
  // Compute ||uh||_Linf(cnc)
  LINF,
  // Compute ||u - uh||_Linf(cnc), with u a user-defined scalar field
  LINF_ERROR,
  // Compute |uh|_H1(cnc)
  SEMI_H1,
  // Compute |u - uh|_H1(cnc), with u a user-defined scalar field
  SEMI_H1_ERROR,
  // Compute ||uh||_H1(cnc)
  H1,
  // Compute ||u - uh||_H1(cnc), with u a user-defined scalar field
  H1_ERROR,
  //         /
  // Compute |    1 dx
  //         / cnc
  AREA,
  //         /
  // Compute |    uh dx
  //         / cnc
  INTEGRAL,
  //         /
  // Compute |    uh dot F dx, with F a user-defined scalar or vector field
  //         / cnc
  DOT_PRODUCT,
  //         /
  // Compute |    f dx, with f a user-defined scalar field at t = 0
  //         / cnc
  INTEGRAL_F
} normType;

class feNorm;

// Create a norm object and perform safety checks.
// This is the recommended way of creating a norm.
// Call within feCheck( ... ) to safely exit if a problem arises.
//
// Although a type is necessary to create a norm, any norm object can call any 
// norm or integral computation routine (just make sure adequate FE spaces are provided).
// The type is used when norms are given to a tim integration scheme, and the generic 
// 'compute' routine is called at each time step.
//
//           type: the norm of integral to compute, see enum above
//         spaces: the FE spaces used to compute the norm. The spaces must be 
//                 defined on the same connectivity (Physical Entity) and be compatible
//                 with the required norm, e.g. an L2 error norm of a scalar field will
//                 be computed using the first FE space only, but the norm of a vector 
//                 field used the first n spaces, with n the space dimension.
//            sol: ptr to the active solution
// scalarSolution: ptr to a scalar field used in the definition of the norm or 
//                 integral to compute, see enum (typically the exact solution)
// vectorSolution: ptr to a vector field used in the definition of the norm or 
//                 integral to compute, see enum (typically the exact solution)
feStatus createNorm(feNorm *&norm,
  normType type,
  const std::vector<feSpace*> &spaces,
  feSolution *sol,
  feFunction *scalarSolution = nullptr,
  feVectorFunction *vectorSolution = nullptr);

//
// Computes a variety of norms and integrals over a geometric connectivity
//
class feNorm
{
protected:
  // The type of norm/integral to compute
  normType _type;
  // FE spaces
  std::vector<feSpace *> _spaces;
  // Convenience ptr to the solution and the geometric connectivity
  const feSolution *_solution;
  const feCncGeo *_cnc;
  // Number of elements on the geometric connectivity
  int _nElm;
  // FE space used to interpolate the geometry
  feSpace *_geoSpace;
  // Quadrature rule of the FE spaces and jacobians
  std::vector<double> _w;
  int _nQuad;
  const std::vector<double> &_J;

  // Addressing vector and local solution for each FE space
  std::vector<std::vector<feInt>> _adr;
  std::vector<std::vector<double>> _localSol;
  // Convenience vectors to evaluate a field at position _pos
  // and physical coordinates of an element
  std::vector<double> _pos;
  std::vector<double> _geoCoord;

  feFunction *_scalarSolution;
  feVectorFunction *_vectorSolution = nullptr;

  bool _cncOnly;

public:
  // Create a norm. To perform safety checks, create a norm using createNorm instead.
  feNorm(normType type, const std::vector<feSpace*> &spaces, feSolution *sol, feFunction *scalarSolution = nullptr,
         feVectorFunction *vectorSolution = nullptr);
  // Create a norm providing only the geometric connectivity and a solution.
  // Used to compute integrals over user-defined field only (not on FE solution).
  feNorm(feCncGeo *cnc, feFunction *scalarSolution = nullptr, feVectorFunction *vectorSolution = nullptr);
  ~feNorm(){};

  // Compute and return the norm or integral matching the type 'type' (see enum above)
  double compute(normType type);
  // Compute and return the norm of integral matching the norm's _type attribute
  double compute();

private:
  void initializeLocalSolutionOnSpace(int iSpace, int iElm);
  double computeLpNorm(int p, bool error = false);
  double computeVectorLpNorm(int p, bool error = false);
  double computeL1Norm(bool error = false);
  double computeL2Norm(bool error = false);
  double computeLInfNorm(bool error = false);
  double computeH1SemiNorm(bool error = false);
  double computeH1Norm(bool error = false);
  double computeArea();
  double computeIntegral();
  double computeIntegralUserFunction();
  double computeIntegralDotProduct();
  feStatus computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                feMesh *mesh, feMetaNumber *refMN,
                                                feSolution *refSol, feMesh *refMesh,
                                                const std::vector<feSpace *> refSpaces,
                                                double &res);
};

#endif
