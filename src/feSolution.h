#ifndef _FESOLUTION_
#define _FESOLUTION_

#include "feMesh.h"
#include "feSpace.h"
#include "feNumber.h"

class feSolution{

protected:
  std::vector<double> _sol;
  std::vector<double> _dsoldt;

  double _c0;
  double _tn; // The current time
  // const std::vector<double>& _coord;

  const std::vector<feSpace*>& _space;
  const std::vector<feSpace*>& _essBC;

  double _t0;
  double _t1;
  int _nTimeSteps;
  double _dt;
public:
	feSolution(feMesh *mesh, const std::vector<feSpace*> &space, const std::vector<feSpace*> &essBC, feMetaNumber *metaNumber)
    : _space(space), _essBC(essBC)//, _coord(mesh->getCoord())
	{
    _sol.resize(metaNumber->getNbDOFs());
    _dsoldt.resize(metaNumber->getNbDOFs());
    for(int i = 0; i < metaNumber->getNbDOFs(); ++i){
      _sol[i] = 0.0;
      _dsoldt[i] = 0.0;
    }

    _c0 = 0.0;
    //_soln 
	};
	~feSolution() {
  }

  void initializeTemporalSolution(double t0, double t1, int nTimeSteps);
  void initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber);
  void initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber);


  int getC0(){ return _c0; }
  double getCurrentTime(){ return _tn; }
  void setCurrentTime(double t){ _tn = t; }

  double getSolAtDOF(int iDOF){ return _sol[iDOF]; }
  double getSolDotAtDOF(int iDOF){ return _dsoldt[iDOF]; }

  void printSol();

};

#endif