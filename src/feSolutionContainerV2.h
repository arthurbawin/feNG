#ifndef _FESOLUTIONCONTAINERV2_
#define _FESOLUTIONCONTAINERV2_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feLinearSystem;

class feSolutionContainerV2
{
protected:
  int _nbInc;
  int _nbDDL;

  double _tC;

  // nSauvegarde
  int _orderBDF;
  int _nStorage; 
  
  std::vector<double> _cBDF;
  std::vector<double> _d;

  std::string _correctionType;
public:
  std::vector<std::vector<double> > _sol;    // U
  std::vector<std::vector<double> > _solDot; // Udot
  std::vector<std::vector<double> > _fResidual; // Udot
  
public:
  feSolutionContainerV2(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber): _orderBDF(orderBDF), _nStorage(nStorage), _tC(tC)
  {
    _nbInc = metaNumber->getNbUnknowns();
    _nbDDL = metaNumber->getNbDOFs();

    _cBDF.resize(_orderBDF+1);

    _d.resize(_nbDDL,0.);

    _sol.resize(_nStorage);
    _solDot.resize(_nStorage);
    _fResidual.resize(_nStorage);

    for (int i=0; i<_nStorage; i++){
      _sol[i].resize(_nbDDL);
      _solDot[i].resize(_nbDDL, 0.);
      _fResidual[i].resize(_nbDDL, 0.);
    }

  };
  
  int getNbDof(){return _nbDDL;};
  int getNbInc(){return _nbInc;};

  std::vector<double> &getCoeffBDF(){return _cBDF;};
  std::vector<double> &getCorrection(){return _d;};

  std::vector<double> &getSol(int n){return _sol[n];};
  std::vector<double> &getSolDot(int n){return _solDot[n];};
  std::vector<double> &getResidual(int n){return _fResidual[n];};

  void copyCurrentSolution(feSolution *currentSolution);
  void copyCurrentSolutionDot(feSolution *currentSolution);

  void computeCoeffBDF(std::vector<double> &k);

  void initializeBDFSolution(feMesh *mesh, feMetaNumber *metaNumber, feSolution *currentSolution);
  void initializeDCSolution(feMesh *mesh, feMetaNumber *metaNumber, feSolution *currentSolution, std::string codeIniDC);
  
  void rotate();
  void invRotate();
  
  void printSolution(int i);


  void computeSolDot(feLinearSystem *linearSystem);

  virtual ~feSolutionContainerV2() {}

};


// ============================================================
// Stationnary
// ============================================================
  //ajouter verif sur nStorage
class feSolutionContainerV2Stat : public feSolutionContainerV2
{
protected:

public:
  feSolutionContainerV2Stat(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber){};
  
  void computeCorrection();

  virtual ~feSolutionContainerV2Stat() {}
  
};


// ============================================================
// BDF
// ============================================================
class feSolutionContainerV2BDF : public feSolutionContainerV2
{
protected:

public:
  feSolutionContainerV2BDF(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber){};

  void computeCorrection();
  
  virtual ~feSolutionContainerV2BDF() {}  
};



// ============================================================
// DC/BDF2
// ============================================================
  //ajouter verification sur orderBDF

class feSolutionContainerV2DC2F : public feSolutionContainerV2
{
protected:
  std::vector<double> _cDC2;

public:
  feSolutionContainerV2DC2F(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _cDC2.resize(2);
  };

  void computeCorrectionCoeff(std::vector<double> &k);

  void computeCorrection(feSolutionContainerV2BDF *solutionContainerBDF1, std::vector<double> &k);//correctionType 
  
  virtual ~feSolutionContainerV2DC2F() {}
  
};


class feSolutionContainerV2DC3F : public feSolutionContainerV2
{
protected:
  std::vector<double> _cDC3;
  std::vector<double> _bDC3;

  std::vector<double> _cStartDC3_t1;
  std::vector<double> _bStartDC3_t1;

public:
  feSolutionContainerV2DC3F(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _cDC3.resize(3);
    _bDC3.resize(3);
  };
  
  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2DC2F *solutionContainerDC2F, std::vector<double> &k);

  void computeStartCorrectionCoeff(std::vector<double> &k);
  void computeStartCorrection(feSolutionContainerV2DC2F *solutionContainerDC2F, std::vector<double> &k);

  virtual ~feSolutionContainerV2DC3F() {}
  
};


class feSolutionContainerV2DC4F : public feSolutionContainerV2
{
protected:
  std::vector<double> _cDC4;
  std::vector<double> _bDC4;
  std::vector<double> _aDC4;

  std::vector<double> _cStartDC4_t1;
  std::vector<double> _bStartDC4_t1;
  std::vector<double> _aStartDC4_t1;

  std::vector<double> _cStartDC4_t2;
  std::vector<double> _bStartDC4_t2;
  std::vector<double> _aStartDC4_t2;

  std::vector<double> _dDC4_t1;
  std::vector<double> _dDC4_t2;

public:
  feSolutionContainerV2DC4F(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _cDC4.resize(4);
    _bDC4.resize(4);
    _aDC4.resize(4);
  };

  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2DC3F *solutionContainerDC3F, std::vector<double> &k);

  void computeStartCorrectionCoeff(std::vector<double> &k);
  void computeStartCorrection(feSolutionContainerV2DC3F *solutionContainerDC3F, std::vector<double> &k);

  void setStartCorrection(int step);
  
  virtual ~feSolutionContainerV2DC4F() {}
  
};

class feSolutionContainerV2DC5F : public feSolutionContainerV2
{
protected:
  std::vector<double> _cDC5;
  std::vector<double> _bDC5;
  std::vector<double> _aDC5;
  std::vector<double> _eDC5;

  std::vector<double> _cStartDC5_t1;
  std::vector<double> _bStartDC5_t1;
  std::vector<double> _aStartDC5_t1;
  std::vector<double> _eStartDC5_t1;

  std::vector<double> _cStartDC5_t2;
  std::vector<double> _bStartDC5_t2;
  std::vector<double> _aStartDC5_t2;
  std::vector<double> _eStartDC5_t2;

  std::vector<double> _cStartDC5_t3;
  std::vector<double> _bStartDC5_t3;
  std::vector<double> _aStartDC5_t3;
  std::vector<double> _eStartDC5_t3;

  std::vector<double> _dDC5_t1;
  std::vector<double> _dDC5_t2;
  std::vector<double> _dDC5_t3;

public:
  feSolutionContainerV2DC5F(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _cDC5.resize(5);
    _bDC5.resize(5);
    _aDC5.resize(5);
    _eDC5.resize(5);
  };

  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2DC4F *solutionContainerDC4F, std::vector<double> &k);
  
  void computeStartCorrectionCoeff(std::vector<double> &k);
  void computeStartCorrection(feSolutionContainerV2DC4F *solutionContainerDC4F, std::vector<double> &k);

  void setStartCorrection(int step);

  virtual ~feSolutionContainerV2DC5F() {}
  
};



// ============================================================
// DC/BDF2
// ============================================================
  //ajouter verification sur orderBDF

class feSolutionContainerV2DC3 : public feSolutionContainerV2
{
protected:
  std::vector<double> _bDC3;

public:
  feSolutionContainerV2DC3(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _bDC3.resize(3);
  };
  
  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2BDF *solutionContainerBDF2, std::vector<double> &k);

  virtual ~feSolutionContainerV2DC3() {}
  
};


class feSolutionContainerV2DC4 : public feSolutionContainerV2
{
protected:
  std::vector<double> _bDC4;
  std::vector<double> _aDC4;

public:
  feSolutionContainerV2DC4(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _bDC4.resize(4);
    _aDC4.resize(4);
  };

  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2DC3 *solutionContainerDC3, std::vector<double> &k);
  
  virtual ~feSolutionContainerV2DC4() {}
  
};

class feSolutionContainerV2DC5 : public feSolutionContainerV2
{
protected:
  std::vector<double> _bDC5;
  std::vector<double> _aDC5;
  std::vector<double> _eDC5;

public:
  feSolutionContainerV2DC5(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _bDC5.resize(5);
    _aDC5.resize(5);
    _eDC5.resize(5);
  };

  void computeCorrectionCoeff(std::vector<double> &k);

  void computeCorrection(feSolutionContainerV2DC4 *solutionContainerDC4, std::vector<double> &k);
  
  virtual ~feSolutionContainerV2DC5() {}
  
};


class feSolutionContainerV2DC6 : public feSolutionContainerV2
{
protected:
  std::vector<double> _bDC6;
  std::vector<double> _aDC6;
  std::vector<double> _eDC6;
  std::vector<double> _fDC6;

public:
  feSolutionContainerV2DC6(int orderBDF, int nStorage, double tC, feMetaNumber *metaNumber, std::string correctionType)
    : feSolutionContainerV2(orderBDF, nStorage, tC, metaNumber)
  {
    _correctionType = correctionType;

    _bDC6.resize(6);
    _aDC6.resize(6);
    _eDC6.resize(6);
    _fDC6.resize(6);
  };

  void computeCorrectionCoeff(std::vector<double> &k);
  void computeCorrection(feSolutionContainerV2DC5 *solutionContainerDC5, std::vector<double> &k);
  
  virtual ~feSolutionContainerV2DC6() {}
  
};



#endif