#ifndef _FETIMEINTEGRATOR_
#define _FETIMEINTEGRATOR_

#include "feMessage.h"
#include "feMesh.h"
#include "feLinearSystem.h"
#include "feBilinearForm.h"
#include "feNumber.h"
// #include "feNorm.h"
#include "feNormV2.h"
#include "feSpace.h"
#include "feExporter.h"
#include "feSolverV2.h"
// #include "feSolutionContainerV2.h"


typedef enum {STATIONARY_V2, BDF1_V2, BDF2_V2, DC2F_V2, DC3F_V2, DC4F_V2, DC5F_V2, DC3_V2, DC4_V2, DC5_V2, DC6_V2} timeIntegratorSchemeV2;



class feTimeIntegrator {
protected:
	feMesh         *_mesh;
	feSolverV2     *_solver;
	feLinearSystem *_linearSystem;
	feMetaNumber   *_metaNumber;
	feSolution     *_currentSolution;

	timeIntegratorSchemeV2 _scheme;
	
	int _nbInc;
	int _nbDDL;

	int    _nTimeStep;
	double _t0;
	double _tEnd;
	double _tC;

	std::vector<double> _dt;
	
	std::string _codeDt;

	std::vector<double> _tG;
	std::vector<double> _t;
	std::vector<double> _k;

	double _gamma;
	double _r;
	double _k0;
	int _nTimeStepM;

	int _nbInterv;
	int _nTimeStepPerInterv;


	int _orderBDF;
	int _nStorage;

	std::vector<feNormV2 *> &_norms;

	feExportData _exportData;


public:
	feTimeIntegrator(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt);

	void initialize();
	void initializeUDOT();

	virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution") = 0;

	virtual std::vector<std::vector<double>> &getNormStat(){};
	virtual std::vector<std::vector<double>> &getNormBDF(){};
  virtual std::vector<std::vector<double>> &getNormDC2(){};
  virtual std::vector<std::vector<double>> &getNormDC3(){};
  virtual std::vector<std::vector<double>> &getNormDC4(){};
  virtual std::vector<std::vector<double>> &getNormDC5(){};
  virtual std::vector<std::vector<double>> &getNormDC6(){};
  virtual std::vector<double> &getINormStat(int i){};
  virtual std::vector<double> &getINormBDF(int i){};
  virtual std::vector<double> &getINormDC2(int i){};
  virtual std::vector<double> &getINormDC3(int i){};
  virtual std::vector<double> &getINormDC4(int i){};
  virtual std::vector<double> &getINormDC5(int i){};
  virtual std::vector<double> &getINormDC6(int i){};

  virtual std::vector<double> &getAllTimes(){return _tG;};
  
  void initializeTime();
	void updateTime(int n);
	void rebootTime();

	void pstClc(feSolutionContainerV2 *solContainer);


	virtual ~feTimeIntegrator() {}
};



// =======================================================
// Stationnary
// =======================================================
class feTimeIntegratorStationary : public feTimeIntegrator
{
protected:
	feSolutionContainerV2Stat *_solutionContainerStat; 
	
	std::vector<std::vector<double>> _resultNormStat;

public:
  feTimeIntegratorStationary(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, double t0, std::vector<feNormV2 *> &norms, feExportData exporter)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, t0, 1, norms, exporter, "Constant")
  {
  	_orderBDF = 0;
		_nStorage = 1;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormStat.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormStat[i].resize(1);
  	}

  	_solutionContainerStat = new feSolutionContainerV2Stat(_orderBDF, _nStorage, _t0, _metaNumber);

  	this->initializeStat();
  }
  
  void initializeStat();

  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormStat(){return _resultNormStat;};
  virtual std::vector<double> &getINormStat(int i){return _resultNormStat[i];};

  virtual ~feTimeIntegratorStationary() {}
};




// =================================================
// BDF
// =================================================
class feTimeIntegratorBDF1 : public feTimeIntegrator
{
protected:
	feSolutionContainerV2BDF *_solutionContainerBDF1; 

	std::vector<std::vector<double>> _resultNormBDF1;

public:
  feTimeIntegratorBDF1(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_orderBDF = 1;
		_nStorage = 2;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF1.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF1[i].resize(_nTimeStep);
  	}

  	_solutionContainerBDF1 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);

  	this->initializeBDF1();
  };
  
  void initializeBDF1();

  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF1;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF1[i];};

  virtual ~feTimeIntegratorBDF1() {}
};




class feTimeIntegratorBDF2 : public feTimeIntegrator
{
protected:
	feSolutionContainerV2BDF *_solutionContainerBDF2; 
	
	std::vector<std::vector<double>> _resultNormBDF2;

public:
  feTimeIntegratorBDF2(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_orderBDF = 2;
		_nStorage = 3;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF2.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF2[i].resize(_nTimeStep);
  	}

  	_solutionContainerBDF2 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);

  	this->initializeBDF2();
  };
  
  void initializeBDF2();

  void startBDF2(std::string startCode);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF2;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF2[i];};
  
  virtual ~feTimeIntegratorBDF2() {}
};



// =================================================
// DC/BDF1
// =================================================
class feTimeIntegratorDC2F : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF1; 
	feSolutionContainerV2DC2F *_solutionContainerDC2F; 

	std::vector<std::vector<double>> _resultNormBDF1;
	std::vector<std::vector<double>> _resultNormDC2F;

public:
  feTimeIntegratorDC2F(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 1;
		_nStorage = 2;
  	
  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF1.resize(_norms.size());
  	_resultNormDC2F.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF1[i].resize(_nTimeStep);
  		_resultNormDC2F[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF1 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC2F = new feSolutionContainerV2DC2F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC2F();
  };

  void initializeDC2F();

  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");
  
  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF1;};
  virtual std::vector<std::vector<double>> &getNormDC2(){return _resultNormDC2F;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF1[i];};
  virtual std::vector<double> &getINormDC2(int i){return _resultNormDC2F[i];};

  virtual ~feTimeIntegratorDC2F() {}
};


class feTimeIntegratorDC3F : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF1; 
	feSolutionContainerV2DC2F *_solutionContainerDC2F; 
	feSolutionContainerV2DC3F *_solutionContainerDC3F;

	std::vector<std::vector<double>> _resultNormBDF1;
	std::vector<std::vector<double>> _resultNormDC2F;
	std::vector<std::vector<double>> _resultNormDC3F;

public:
  feTimeIntegratorDC3F(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 1;
		_nStorage = 3;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);
  	
  	_resultNormBDF1.resize(_norms.size());
  	_resultNormDC2F.resize(_norms.size());
  	_resultNormDC3F.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF1[i].resize(_nTimeStep);
  		_resultNormDC2F[i].resize(_nTimeStep);
			_resultNormDC3F[i].resize(_nTimeStep);
  	}

  	_solutionContainerBDF1 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC2F = new feSolutionContainerV2DC2F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC3F = new feSolutionContainerV2DC3F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC3F();
  };

  void initializeDC3F();

  void startDC3F(std::string startCode, std::string codeIniDC);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF1;};
  virtual std::vector<std::vector<double>> &getNormDC2(){return _resultNormDC2F;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3F;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF1[i];};
  virtual std::vector<double> &getINormDC2(int i){return _resultNormDC2F[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3F[i];};

  virtual ~feTimeIntegratorDC3F() {}
};


class feTimeIntegratorDC4F : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF  *_solutionContainerBDF1; 
	feSolutionContainerV2DC2F *_solutionContainerDC2F; 
	feSolutionContainerV2DC3F *_solutionContainerDC3F;
	feSolutionContainerV2DC4F *_solutionContainerDC4F;

	std::vector<std::vector<double>> _resultNormBDF1;
	std::vector<std::vector<double>> _resultNormDC2F;
	std::vector<std::vector<double>> _resultNormDC3F;
	std::vector<std::vector<double>> _resultNormDC4F;

public:
  feTimeIntegratorDC4F(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 1;
		_nStorage = 4;
  	
  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF1.resize(_norms.size());
  	_resultNormDC2F.resize(_norms.size());
  	_resultNormDC3F.resize(_norms.size());
  	_resultNormDC4F.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF1[i].resize(_nTimeStep);
  		_resultNormDC2F[i].resize(_nTimeStep);
			_resultNormDC3F[i].resize(_nTimeStep);
			_resultNormDC4F[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF1 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC2F = new feSolutionContainerV2DC2F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC3F = new feSolutionContainerV2DC3F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC4F = new feSolutionContainerV2DC4F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC4F();
  };

  void initializeDC4F();

  void startDC4F(std::string startCode, std::string codeIniDC);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF1;};
  virtual std::vector<std::vector<double>> &getNormDC2(){return _resultNormDC2F;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3F;};
  virtual std::vector<std::vector<double>> &getNormDC4(){return _resultNormDC4F;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF1[i];};
  virtual std::vector<double> &getINormDC2(int i){return _resultNormDC2F[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3F[i];};
  virtual std::vector<double> &getINormDC4(int i){return _resultNormDC4F[i];};
  
  virtual ~feTimeIntegratorDC4F() {}
};

class feTimeIntegratorDC5F : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF1; 
	feSolutionContainerV2DC2F *_solutionContainerDC2F; 
	feSolutionContainerV2DC3F *_solutionContainerDC3F;
	feSolutionContainerV2DC4F *_solutionContainerDC4F;
	feSolutionContainerV2DC5F *_solutionContainerDC5F;

	std::vector<std::vector<double>> _resultNormBDF1;
	std::vector<std::vector<double>> _resultNormDC2F;
	std::vector<std::vector<double>> _resultNormDC3F;
	std::vector<std::vector<double>> _resultNormDC4F;
	std::vector<std::vector<double>> _resultNormDC5F;

public:
  feTimeIntegratorDC5F(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
  		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 1;
		_nStorage = 5;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF1.resize(_norms.size());
  	_resultNormDC2F.resize(_norms.size());
  	_resultNormDC3F.resize(_norms.size());
  	_resultNormDC4F.resize(_norms.size());
  	_resultNormDC5F.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF1[i].resize(_nTimeStep);
  		_resultNormDC2F[i].resize(_nTimeStep);
			_resultNormDC3F[i].resize(_nTimeStep);
			_resultNormDC4F[i].resize(_nTimeStep);
			_resultNormDC5F[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF1 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC2F = new feSolutionContainerV2DC2F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC3F = new feSolutionContainerV2DC3F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC4F = new feSolutionContainerV2DC4F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC5F = new feSolutionContainerV2DC5F(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC5F();
  };

  void initializeDC5F();

  void startDC5F(std::string startCode, std::string codeIniDC);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF1;};
  virtual std::vector<std::vector<double>> &getNormDC2(){return _resultNormDC2F;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3F;};
  virtual std::vector<std::vector<double>> &getNormDC4(){return _resultNormDC4F;};
  virtual std::vector<std::vector<double>> &getNormDC5(){return _resultNormDC5F;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF1[i];};
  virtual std::vector<double> &getINormDC2(int i){return _resultNormDC2F[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3F[i];};
  virtual std::vector<double> &getINormDC4(int i){return _resultNormDC4F[i];};
  virtual std::vector<double> &getINormDC5(int i){return _resultNormDC5F[i];};

  virtual ~feTimeIntegratorDC5F() {}
};


// ==================================================
// DC/BDF2
// ==================================================
class feTimeIntegratorDC3 : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF2; 
	feSolutionContainerV2DC3 *_solutionContainerDC3;

	std::vector<std::vector<double>> _resultNormBDF2;
	std::vector<std::vector<double>> _resultNormDC3;

public:
  feTimeIntegratorDC3(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
   		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 2;
		_nStorage = 3;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);
  	
  	_resultNormBDF2.resize(_norms.size());
  	_resultNormDC3.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF2[i].resize(_nTimeStep);
			_resultNormDC3[i].resize(_nTimeStep);
  	}

  	_solutionContainerBDF2 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC3  = new feSolutionContainerV2DC3(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC3();
  };

  void initializeDC3();

  void startDC3(std::string startCode);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF2;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF2[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3[i];};
  
  virtual ~feTimeIntegratorDC3() {}
};



class feTimeIntegratorDC4 : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF2; 
	feSolutionContainerV2DC3 *_solutionContainerDC3;
	feSolutionContainerV2DC4 *_solutionContainerDC4;

	std::vector<std::vector<double>> _resultNormBDF2;
	std::vector<std::vector<double>> _resultNormDC3;
	std::vector<std::vector<double>> _resultNormDC4;

public:
  feTimeIntegratorDC4(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
   		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 2;	
		_nStorage = 4;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF2.resize(_norms.size());
  	_resultNormDC3.resize(_norms.size());
  	_resultNormDC4.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF2[i].resize(_nTimeStep);
			_resultNormDC3[i].resize(_nTimeStep);
			_resultNormDC4[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF2 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC3  = new feSolutionContainerV2DC3(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC4  = new feSolutionContainerV2DC4(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC4();
  };
  
  void initializeDC4();

  void startDC4(std::string startCode);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF2;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3;};
  virtual std::vector<std::vector<double>> &getNormDC4(){return _resultNormDC4;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF2[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3[i];};
  virtual std::vector<double> &getINormDC4(int i){return _resultNormDC4[i];};

  virtual ~feTimeIntegratorDC4() {}
};



class feTimeIntegratorDC5 : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF2; 
	feSolutionContainerV2DC3 *_solutionContainerDC3;
	feSolutionContainerV2DC4 *_solutionContainerDC4;
	feSolutionContainerV2DC5 *_solutionContainerDC5;

	std::vector<std::vector<double>> _resultNormBDF2;
	std::vector<std::vector<double>> _resultNormDC3;
	std::vector<std::vector<double>> _resultNormDC4;
	std::vector<std::vector<double>> _resultNormDC5;

public:
  feTimeIntegratorDC5(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
   		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 2;
		_nStorage = 5;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF2.resize(_norms.size());
  	_resultNormDC3.resize(_norms.size());
  	_resultNormDC4.resize(_norms.size());
  	_resultNormDC5.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF2[i].resize(_nTimeStep);
			_resultNormDC3[i].resize(_nTimeStep);
			_resultNormDC4[i].resize(_nTimeStep);
			_resultNormDC5[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF2 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC3  = new feSolutionContainerV2DC3(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC4  = new feSolutionContainerV2DC4(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC5  = new feSolutionContainerV2DC5(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC5();
  };
  
  void initializeDC5();

  void startDC5(std::string startCode);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF2;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3;};
  virtual std::vector<std::vector<double>> &getNormDC4(){return _resultNormDC4;};
  virtual std::vector<std::vector<double>> &getNormDC5(){return _resultNormDC5;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF2[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3[i];};
  virtual std::vector<double> &getINormDC4(int i){return _resultNormDC4[i];};
  virtual std::vector<double> &getINormDC5(int i){return _resultNormDC5[i];};

  virtual ~feTimeIntegratorDC5() {}
};


class feTimeIntegratorDC6 : public feTimeIntegrator
{
protected:
	std::string _correctionType;

	feSolutionContainerV2BDF *_solutionContainerBDF2; 
	feSolutionContainerV2DC3 *_solutionContainerDC3;
	feSolutionContainerV2DC4 *_solutionContainerDC4;
	feSolutionContainerV2DC5 *_solutionContainerDC5;
	feSolutionContainerV2DC6 *_solutionContainerDC6;

	std::vector<std::vector<double>> _resultNormBDF2;
	std::vector<std::vector<double>> _resultNormDC3;
	std::vector<std::vector<double>> _resultNormDC4;
	std::vector<std::vector<double>> _resultNormDC5;
	std::vector<std::vector<double>> _resultNormDC6;

public:
  feTimeIntegratorDC6(feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt)
   		:feTimeIntegrator(mesh, solver, linearSystem, metaNumber, currentSolution, scheme, t0, tEnd, nTimeStep, norms, exporter, codeDt)
  {
  	_correctionType = correcType;

  	_orderBDF = 2;
		_nStorage = 6;

  	_k.resize(_nStorage,0.0);
  	_t.resize(_nStorage,0.0);

  	_resultNormBDF2.resize(_norms.size());
  	_resultNormDC3.resize(_norms.size());
  	_resultNormDC4.resize(_norms.size());
  	_resultNormDC5.resize(_norms.size());
  	_resultNormDC6.resize(_norms.size());

  	for (size_t i=0; i<_norms.size();i++){
  		_resultNormBDF2[i].resize(_nTimeStep);
			_resultNormDC3[i].resize(_nTimeStep);
			_resultNormDC4[i].resize(_nTimeStep);
			_resultNormDC5[i].resize(_nTimeStep);
			_resultNormDC6[i].resize(_nTimeStep);
  	}
  	
  	_solutionContainerBDF2 = new feSolutionContainerV2BDF(_orderBDF, _nStorage, _t0, _metaNumber);
  	_solutionContainerDC3  = new feSolutionContainerV2DC3(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC4  = new feSolutionContainerV2DC4(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC5  = new feSolutionContainerV2DC5(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);
  	_solutionContainerDC6  = new feSolutionContainerV2DC6(_orderBDF, _nStorage, _t0, _metaNumber, _correctionType);

  	this->initializeDC6();
  };
  
  void initializeDC6();

  void startDC6(std::string startCode);
  virtual feStatus makeSteps(std::string startCode="Analytique", std::string codeIniDC="fromPreviousSolution");

  virtual std::vector<std::vector<double>> &getNormBDF(){return _resultNormBDF2;};
  virtual std::vector<std::vector<double>> &getNormDC3(){return _resultNormDC3;};
  virtual std::vector<std::vector<double>> &getNormDC4(){return _resultNormDC4;};
  virtual std::vector<std::vector<double>> &getNormDC5(){return _resultNormDC5;};
  virtual std::vector<std::vector<double>> &getNormDC6(){return _resultNormDC6;};
  virtual std::vector<double> &getINormBDF(int i){return _resultNormBDF2[i];};
  virtual std::vector<double> &getINormDC3(int i){return _resultNormDC3[i];};
  virtual std::vector<double> &getINormDC4(int i){return _resultNormDC4[i];};
  virtual std::vector<double> &getINormDC5(int i){return _resultNormDC5[i];};
  virtual std::vector<double> &getINormDC6(int i){return _resultNormDC6[i];};

  virtual ~feTimeIntegratorDC6() {}
};

feStatus createTimeIntegratorV2(feTimeIntegrator *&timeIntegrator, feMesh *mesh, feSolverV2 *solver, feLinearSystem *linearSystem, feMetaNumber *metaNumber, feSolution *currentSolution, timeIntegratorSchemeV2 scheme, std::string correcType, double t0, double tEnd, int nTimeStep, std::vector<feNormV2 *> &norms, feExportData exporter, std::string codeDt);

#endif