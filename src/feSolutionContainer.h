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
  std::vector<double> _d;
public:
  std::vector<double> _cn;
  std::vector<std::vector<double>> _sol; // U
  std::vector<std::vector<double>> _fResidual; // F
public:
	feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber);
	virtual ~feSolutionContainer() {}

  std::vector<double> &getTime(){ return _t; }

  void initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber);
  void rotate(double dt);
  void setSol(int iSol, std::vector<double> sol){ _sol[iSol] = sol; }
  // double getSol(int iSol, int iDOF){ return _sol[iSol][iDOF]; }
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

void initializeBDF2(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, 
  feSolutionBDF2 *solBDF2);
void initializeDC3F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, 
  feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3);
void initializeDC4F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, 
  feSolutionBDF2 *solBDF2, feSolutionDCF *solDC3, feSolutionDCF *solDC4);

#endif