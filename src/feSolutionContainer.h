#ifndef _FESOLUTIONCONTAINER_
#define _FESOLUTIONCONTAINER_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feLinearSystem;

class feSolutionContainer {
protected:
  int _nDofs; // NDDL
  int _nSol;
  std::vector<double> _t; // TEMPS
public:
  std::vector<double> _cn;
  std::vector<double> _d;
  std::vector<std::vector<double> > _sol; // U
  std::vector<std::vector<double> > _fResidual; // F
public:
  feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber);
  virtual ~feSolutionContainer() {}

  int getNbDOFs() { return _nDofs; }
  void setNbDOFs(int nDofs) { _nDofs = nDofs; }

  int getNbSol() { return _nSol; }
  std::vector<double> &getTime() { return _t; }
  std::vector<double> &getSolution(int iSol) { return _sol[iSol]; }

  void initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber);
  void rotate(double dt);
  void setSol(int iSol, std::vector<double> sol) { _sol[iSol] = sol; }
  void setSolAtDOF(int iSol, int iDOF, double val) { _sol[iSol][iDOF] = val; }
  // double getSol(int iSol, int iDOF){ return _sol[iSol][iDOF]; }
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem){};
};

class feStationarySolution : public feSolutionContainer {
protected:
public:
  feStationarySolution(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feStationarySolution() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF1 : public feSolutionContainer {
protected:
public:
  feSolutionBDF1(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionBDF1() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF2 : public feSolutionContainer {
protected:
public:
  feSolutionBDF2(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionBDF2() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionDCF : public feSolutionContainer {
protected:
public:
  feSolutionDCF(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionDCF() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

void initializeBDF2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2);
void initializeBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1);
void initializeBDF2withBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                            feSolutionBDF1 *solBDF1, feSolutionBDF2 *solBDF2);
void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3);
void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4);

#endif