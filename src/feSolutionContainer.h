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
  double getSol(int iSol, int iDOF) { return _sol[iSol][iDOF]; }
  double getRes(int iSol, int iDOF) { return _fResidual[iSol][iDOF]; }
  std::vector<double> &getResidual(int iSol) { return _fResidual[iSol]; }
  void setResidual(int iSol, std::vector<double> res) { _fResidual[iSol] = res; }
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem){};
};

class feStationarySolution : public feSolutionContainer
{
protected:
public:
  feStationarySolution(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feStationarySolution() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

//////////////////////////////////////////////////////////
// Container with solutions u^n , u^n-1  but no correction 
class feSolutionBDF1 : public feSolutionContainer
{
protected:
public:
  feSolutionBDF1(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionBDF1() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};
//////////////////////////////////////////////////////////
// Container with solutions u^n , u^n-1, u^n-2  but no correction 
class feSolutionBDF2 : public feSolutionContainer
{
protected:
public:
  feSolutionBDF2(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionBDF2() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

//////////////////////////////////////////////////////////
// Container with solutions u^n , u^n-1, u^n-2  and a correction _d
//used for DC3, DC4, DC5
class feSolutionDCF : public feSolutionContainer
{
protected:
public:
  feSolutionDCF(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionDCF() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

//////////////////////////////////////////////////////////
// Container with solutions u^n , u^n-1 and a correction _d
// Used for DC2F, DC3F, DC4F, DC5F
class feSolutionDC2F : public feSolutionContainer     
{
protected:
public:
  feSolutionDC2F(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  virtual ~feSolutionDC2F() {}
  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
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
                             feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3);
void initializeDC3F_centered(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC2F, feSolutionDC2F *solDC3);
void initializeDC4(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3, feSolutionDCF *solDC4);
void initializeDC5(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh, feSolutionBDF2 *solBDF2,
                   feSolutionDCF *solDC3, feSolutionDCF *solDC4, feSolutionDCF *solDC5);
void initializeDC5F(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4F, feSolutionDC2F *solDC5);
/////////////////////////////////////////////////////////////////////////////////////
// Set the correction for the first value of DC4F u^1, u^2, u^3 of order 4 .
void initializeDC4F_begining(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4, std::string time );
/////////////////////////////////////////////////////////////////////////////////////
// Set the correction for the first value of DC5F u^1, u^2, u^3, u^4 of order 5 .
void initializeDC5F_begining(feSolution *sol, feMetaNumber *metaNumber, feMesh *mesh,
                             feSolutionDC2F *solDC3F, feSolutionDC2F *solDC4F, feSolutionDC2F *solDC5, std::string time );


#endif