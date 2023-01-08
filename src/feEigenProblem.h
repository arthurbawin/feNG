#ifndef _FEEIGENPROBLEM_
#define _FEEIGENPROBLEM_

#include "feMessage.h"
#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"
#include "feSolutionContainer.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

#if defined(HAVE_SLEPC)
#include "slepceps.h"
#endif

/* Supported eigen solvers */
typedef enum { SLEPC } eigenSolverType;

typedef struct eigenPairType{
  double valReal;
  double valImag;
#if defined(HAVE_PETSC)
  Vec vecReal;
  Vec vecImag;
#endif
} eigenPair;

class feEigenProblem
{
protected:
  size_t _numAForms; // The number of bilinear forms to assemble matrix A in A*x = lambda*B*x
  size_t _numBForms; // The number of bilinear forms to assemble matrix B in A*x = lambda*B*x
  std::vector<feBilinearForm *> _AForms;
  std::vector<feBilinearForm *> _BForms;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  int _argc;
  char **_argv;

#if defined(HAVE_SLEPC)
  PetscInt _nInc;
  PetscInt _nDofs;
  Mat _A,_B;
  EPS _eps; // Eigenproblem solver context
  ST _st;
  KSP _ksp;
  EPSType _type;
  PetscReal _tol;
#endif

  std::vector<eigenPair> _eigenPairs;

public:
  feEigenProblem(int argc, char **argv,
                 std::vector<feBilinearForm *> lhsForms,
                 std::vector<feBilinearForm *> rhsForms,
                 feMetaNumber *metaNumber,
                 feMesh *mesh);

  ~feEigenProblem();

  // void initialize(){};
  void setMatricesToZero();
  void assembleLHSMatrix(feSolution *sol);
  void assembleRHSMatrix(feSolution *sol);
  void assemble(feSolution *sol);
  void solve(feSolution *sol);
  void viewLHSMatrix();
  void viewRHSMatrix();
  std::vector<feBilinearForm *> getLHSForms() { return _AForms; };
  std::vector<feBilinearForm *> getRHSForms() { return _BForms; };

  size_t getNumConvergedPairs() { return _eigenPairs.size(); };
  eigenPair getEigenpair(size_t i) { return _eigenPairs[i]; };

  void setEigenmodeAsActiveSolution(feSolution *sol, size_t i, std::string fieldID);
};

void slepcInitialize(int argc, char **argv);
void slepcFinalize();

feStatus createEigenProblem(feEigenProblem *&eigenProblem,
                            eigenSolverType type,
                            std::vector<feSpace *> allFESpaces,
                            std::vector<feBilinearForm *> lhsForms,
                            std::vector<feBilinearForm *> rhsForms,
                            feMetaNumber *metaNumber,
                            feMesh *mesh, int argc = 0, char **argv = nullptr);

#endif