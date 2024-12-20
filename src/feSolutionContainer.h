#ifndef _FESOLUTIONCONTAINER_
#define _FESOLUTIONCONTAINER_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feLinearSystem;

class feSolutionContainer
{
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
  feSolutionContainer(){};
  feSolutionContainer(int nSol, double tn, int nDOF);
  feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber);
  virtual ~feSolutionContainer() {}

  void copy(const feSolutionContainer &other);

  int getNbDOFs() { return _nDofs; }
  void setNbDOFs(int nDofs) { _nDofs = nDofs; }

  int getNbSol() { return _nSol; }
  std::vector<double> &getTime() { return _t; }
  std::vector<double> &getSolution(int iSol) { return _sol[iSol]; }

  void initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber);
  void rotate(double dt);
  void setSol(int iSol, std::vector<double> sol) { _sol[iSol] = sol; }
  void setSolAtDOF(int iSol, int iDOF, double val) { _sol[iSol][iDOF] = val; }
  double getSol(int iSol, int iDOF) { return _sol[iSol][iDOF]; }
  double getRes(int iSol, int iDOF) { return _fResidual[iSol][iDOF]; }
  std::vector<double> &getResidual(int iSol) { return _fResidual[iSol]; }
  void setResidual(int iSol, std::vector<double> res) { _fResidual[iSol] = res; }
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feStationarySolution : public feSolutionContainer
{
public:
  feStationarySolution(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  ~feStationarySolution() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF1 : public feSolutionContainer
{
public:
  feSolutionBDF1(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  ~feSolutionBDF1() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF2 : public feSolutionContainer
{
public:
  feSolutionBDF2(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  ~feSolutionBDF2() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionDCF : public feSolutionContainer
{
public:
  feSolutionDCF(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  ~feSolutionDCF() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionDC2F : public feSolutionContainer
{
public:
  feSolutionDC2F(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  ~feSolutionDC2F() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

void initializeBDF2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2);
void initializeBDF1(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1);
void initializeDC2F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF1 *solBDF1, feSolutionDC2F *solDC2F);
void initializeDC3(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3);
void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3);
void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                    feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4);
void initializeDC3F_centered(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3);

#endif