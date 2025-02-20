#ifndef _FENORM_
#define _FENORM_

#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feNewRecovery.h"

typedef enum {
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
  L2_ERROR_ESTIMATE,
  VECTOR_L2_ERROR,
  // Compute ||uh||_Linf(cnc)
  LINF,
  // Compute ||u - uh||_Linf(cnc), with u a user-defined scalar field
  LINF_ERROR,
  VECTOR_LINF_ERROR,
  // Compute |uh|_H1(cnc)
  SEMI_H1,
  // Compute |u - uh|_H1(cnc), with u a user-defined scalar field
  SEMI_H1_ERROR,
  SEMI_H1_ERROR_ESTIMATE,
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
  INTEGRAL_F,
  PRESSURE_LIFT_FORCE,
  PRESSURE_DRAG_FORCE,
  VISCOUS_LIFT_FORCE,
  VISCOUS_DRAG_FORCE
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
feStatus createNorm(feNorm *&norm, normType type, const std::vector<feSpace *> &spaces,
                    feSolution *sol, feFunction *scalarSolution = nullptr,
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

  // A reconstruction of the solution (and derivatives) used to estimate the error
  feNewRecovery *_rec = nullptr;

  // Addressing vector and local solution for each FE space
  std::vector<std::vector<feInt> > _adr;
  std::vector<std::vector<double> > _localSol;
  // Convenience vectors to evaluate a field at position _pos
  // and physical coordinates of an element
  std::vector<double> _pos;
  std::vector<double> _geoCoord;

  feFunction *_scalarSolution;
  feVectorFunction *_vectorSolution = nullptr;

  bool _cncOnly;
  // Plot error in Gmsh .pos file?
  bool _plotErrorToFile = false;
  std::string _errorPlotFileName = "errorOnElements.pos";

  // To ignore the top fraction of elements error
  bool _ignoreHighestError = false;
  double _fractionToIgnore = 0.;

  double _viscosity = 1.0;

public:
  // Create a norm. To perform safety checks, create a norm using createNorm instead.
  feNorm(normType type, const std::vector<feSpace *> &spaces, feSolution *sol,
         feFunction *scalarSolution = nullptr, feVectorFunction *vectorSolution = nullptr);
  // Create a norm providing only the geometric connectivity and a solution.
  // Used to compute integrals over user-defined field only (not on FE solution).
  feNorm(feCncGeo *cnc, feFunction *scalarSolution = nullptr,
         feVectorFunction *vectorSolution = nullptr);
  ~feNorm(){};

  // Assign a reconstructed solution
  void setRecovery(feNewRecovery *recovery){ _rec = recovery; };
  feNewRecovery *getRecovery(){ return _rec; };

  void setErrorPlotFlag(bool flag) { _plotErrorToFile = flag; };
  void setErrorPlotFilename(std::string name) { _errorPlotFileName = name; };

  void setFractionToIgnore(double fraction)
  {
    _ignoreHighestError = true;
    _fractionToIgnore = fmin(fmax(fraction, 0.), 1.);
  }

  void setViscosity(double viscosity){ _viscosity = viscosity; };

  // Compute and return the norm or integral matching the type 'type' (see enum above)
  double compute(normType type);
  // Compute and return the norm of integral matching the norm's _type attribute
  double compute();

  double computeSquaredErrorOnElement(int iElm);
  double computeSquaredErrorFromEstimatorOnElement(int iElm, bool useAverageEvaluations);
  double computeLpNormOnElement(int p, bool error, int iElm);
  double computeLpErrorExactVsEstimator(int p);
  double computeLInfErrorExactVsEstimator();
  double computeLInfErrorExactVsEstimatorAtVertices();

  // Bounding box is given as [xmin, xmax, ymin, ymax]
  double computeSemiH1ErrorExactVsEstimator(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeSemiH1ErrorExactVsEstimator_Linf(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorHessianExactVsEstimator(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorHessianExactVsEstimator_Linf(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorThirdDerivativesExactVsEstimator(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorThirdDerivativesExactVsEstimator_Linf(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorFourthDerivativesExactVsEstimator(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);
  double computeErrorFourthDerivativesExactVsEstimator_Linf(bool excludeBoundary, double boundaryWidth, double boundingBox[4]);

  // iComponent = 0 : e_x
  // iComponent = 1 : e_y
  double computeForcesFromLagrangeMultiplier(int iComponent = 0);
  double computeViscousDrag(double viscosity, feNewRecovery *recU, feNewRecovery *recV);

  double computeLpErrorFromTransferredSolution(int p, feSolution *otherSol);

  // void computeInterpolationErrorGradientRochery(const int whichElements[2],
  //                                               const int whichControlPoint_localTag[2],
  //                                               double gradient[2]);

  // Implemented in feMinimizeInterpolation.cpp
  void computeInterpolationErrorGradient(const std::vector<int> &whichElements,
                                         const std::vector<int> &whichControlPoint_localTag,
                                         double gradient[2]);

private:
  void initializeLocalSolutionOnSpace(int iSpace, int iElm);
  double computeLpNorm(int p, bool error = false);
  double computeVectorLpNorm(int p, bool error = false);
  double computeLpErrorEstimator(int p);
  double computeL1Norm(bool error = false);
  double computeL2Norm(bool error = false);
  double computeLInfNorm(bool error = false);
  double computeVectorLInfNorm(bool error = false);
  double computeH1SemiNorm(bool error = false);
  double computeH1SemiNormErrorEstimator();
  double computeH1Norm(bool error = false);
  double computeArea();
  double computeIntegral();
  double computeIntegralUserFunction();
  double computeIntegralDotProduct();
  double computePressureLift();
  double computePressureDrag();
  double computeViscousLift();
  feStatus computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                feMesh *mesh, feMetaNumber *refMN,
                                                feSolution *refSol, feMesh *refMesh,
                                                const std::vector<feSpace *> refSpaces,
                                                double &res);
};

#endif
