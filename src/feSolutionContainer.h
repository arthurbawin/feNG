#ifndef _FESOLUTIONCONTAINER_
#define _FESOLUTIONCONTAINER_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feLinearSystem;

//
// A solution container is a collection of N solution arrays,
// as well as their time derivatives, used to perform time integration.
//
// A stationary solver uses a trivial container with a single solution,
// but more involved time integrators (such as "deferred correction" (DC) schemes)
// require combining multiple solutions of different accuracy.
//
class feSolutionContainer
{
protected:
  int _nDofs;
  int _nSol;
  std::vector<double> _t;
  std::vector<double> _deltaT;

  // TODO: Set these to protected
public:
  std::vector<double> _cn; // BDF coefficients, to move to timeintegrator?
  std::vector<double> _d;  // DC correction or residual, to keep

  std::vector<std::vector<double>> _sol; // U
  std::vector<std::vector<double>> _solDot; // Udot
  std::vector<std::vector<double>> _fResidual; // F

public:
  feSolutionContainer(){};
  feSolutionContainer(int nSol, double tn, int nDOF);
  feSolutionContainer(int nSol, double tn, feMetaNumber *metaNumber);
  virtual ~feSolutionContainer(){};

  int getNbDOFs() const { return _nDofs; }
  void setNbDOFs(int nDofs) { _nDofs = nDofs; }
  int getNbSol() const { return _nSol; }
  double getC0() const { return _cn[0]; }
  const std::vector<double> &getTime() const { return _t; }
  const std::vector<double> &getTimeDifferences() const { return _deltaT; }
  const std::vector<double> &getSolution(int i) { return _sol[i]; }
  const std::vector<double> &getSolution(int i) const { return _sol[i]; }
  const std::vector<double> &getSolutionDot(int i) const { return _solDot[i]; }
  const std::vector<double> &getFResidual(int i) const { return _fResidual[i]; }

  void setTime(const int i, const double time) { _t[i] = time; }
  void setTimeDifference(const int i, const double dt) { _deltaT[i] = dt; }

  void setCurrentSolution(const feSolution *other);
  void setCurrentSolutionDot(const feSolution *other);

  void setSolution(const std::vector<double> &other, const int i) { _sol[i] = other; }
  void setSolutionDot(const std::vector<double> &other, const int i) { _solDot[i] = other; }
  void setResidual(const std::vector<double> &other, const int i) { _fResidual[i] = other; }

  // Compute coefficients for BDF methods of order 0, 1, 2 (0 is steady-state)
  void computeBDFCoefficients(const int order, const std::vector<double> &deltaT);
  const std::vector<double> &getBDFCoefficients() const { return _cn; }
  void setBDFCoefficients(const std::vector<double> &coeff);

  // Compute time derivative of current solution using BDF coefficients
  // and system residual if using DC scheme
  void computeCurrentSolDot(feLinearSystem *linearSystem);

  void rotate(const double dt);
  void invRotate(const double dt);
  void rotateWithoutTime();

  virtual void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);

  /////////////
  void NaNify();
  /////////////

  ///////////////////////////////////////////////////////////////////////////
  // Everything below should be reworked

  void initialize(feSolution *sol, feMesh *mesh);
  
  void setSol(int iSol, std::vector<double> sol) { _sol[iSol] = sol; }
  void setSolAtDOF(int iSol, int iDOF, double val) { _sol[iSol][iDOF] = val; }

  double getSol(int iSol, int iDOF) const { return _sol[iSol][iDOF]; }
  double getRes(int iSol, int iDOF) const { return _fResidual[iSol][iDOF]; }
  std::vector<double> &getResidual(int iSol) { return _fResidual[iSol]; }
  void setResidual(int iSol, std::vector<double> res) { _fResidual[iSol] = res; }
};

// Currently, the only specialization of the derived class BDFContainer
// is the computation of the time derivative with precomputed coefficients.
// But the BDF coefficients are stored in the base class...
// This should be reworked.
class BDFContainer : public feSolutionContainer
{
public:
  BDFContainer(int nSol, double tn, int nDOF)
    : feSolutionContainer(nSol, tn, nDOF){};

  // Copy assignment from base class
  BDFContainer& operator=(const feSolutionContainer &other)
  {
    feSolutionContainer::operator=(other);
    return *this;
  }
  
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class DCBDFContainer : public BDFContainer
{
public:
  // The order of the associated DC method
  int _order;
  // The order of the BDF method it is based on (1 or 2)
  int _orderBDF;

  // Coefficients of the DC correction, i.e.,
  // divided differences coefficients.
  // These are the a_ij, b_ij, c_ij, etc. coefficients
  // in Baptiste Berlioux's master thesis.
  //
  // There are _order - 1 sets of coefficients for BDF1-based methods
  // and _order - 2 sets of coefficients for BDF2-based methods
  std::vector<std::vector<double>> _correctionCoefficients;

public:
  DCBDFContainer(const int order,
                 const int orderBDF,
                 const int nSol,
                 const double tn,
                 const int nDOF)
    : BDFContainer(nSol, tn, nDOF), _order(order), _orderBDF(orderBDF)
  {
    if(_orderBDF == 1)
      _correctionCoefficients.resize(_order-1);
    if(_orderBDF == 2)
      _correctionCoefficients.resize(_order-2);

    for(auto &vec : _correctionCoefficients)
      vec.resize(_order);
  };

  // Copy assignment from base class
  DCBDFContainer& operator=(const feSolutionContainer &other)
  {
    BDFContainer::operator=(other);
    return *this;
  }
  
  void computeDCCoefficients(const std::vector<double> &deltaT);
  void computeFirstDCCoefficients(const std::vector<double> &deltaT);
  void computeDCCorrection(const BDFContainer &bdfcontainer,
                           const std::vector<double> &deltaT);
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

//==============================================
// Below is deprecated
//==============================================
class feStationarySolution : public feSolutionContainer
{
public:
  feStationarySolution(int nSol, double tn, int nDOF)
    : feSolutionContainer(nSol, tn, nDOF){};
  feStationarySolution(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  // ~feStationarySolution() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF1 : public feSolutionContainer
{
public:
  feSolutionBDF1(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  // ~feSolutionBDF1() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionBDF2 : public feSolutionContainer
{
public:
  feSolutionBDF2(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  // ~feSolutionBDF2() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionDCF : public feSolutionContainer
{
public:
  feSolutionDCF(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  // ~feSolutionDCF() {}
  void computeSolTimeDerivative(feSolution *sol, feLinearSystem *linearSystem);
};

class feSolutionDC2F : public feSolutionContainer
{
public:
  feSolutionDC2F(int nSol, double tn, feMetaNumber *metaNumber)
    : feSolutionContainer(nSol, tn, metaNumber){};
  // ~feSolutionDC2F() {}
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