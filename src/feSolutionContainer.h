#ifndef _FESOLUTIONCONTAINER_
#define _FESOLUTIONCONTAINER_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feSolutionContainer{
protected:
  int _nDofs; // NDDL
  int _nSol;
  std::vector<double> _t; // TEMPS
  std::vector<std::vector<double>> _sol; // U
  std::vector<std::vector<double>> _fResidual; // F
  std::vector<double> _cn;
  std::vector<double> _d;
public:
	feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber);
	virtual ~feSolutionContainer() {}

  void initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber);
  void setSol(int iSol, std::vector<double> sol){ _sol[iSol] = sol; }
  virtual void computeSolTimeDerivative(feSolution *sol){};
};

class feStationarySolution : public feSolutionContainer{
protected:
public:
  feStationarySolution(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol,tn,metaNumber){};
  virtual ~feStationarySolution() {}
  virtual void computeSolTimeDerivative(feSolution *sol);
};

class feSolutionBDF2 : public feSolutionContainer{
protected:
public:
  feSolutionBDF2(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol,tn,metaNumber){};
  virtual ~feSolutionBDF2() {}
  virtual void computeSolTimeDerivative(feSolution *sol);
};

class feSolutionDCF : public feSolutionContainer{
protected:
public:
  feSolutionDCF(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol,tn,metaNumber){};
  virtual ~feSolutionDCF() {}
  virtual void computeSolTimeDerivative(feSolution *sol);
};

#endif