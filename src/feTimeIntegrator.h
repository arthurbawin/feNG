#ifndef _FETIMEINTEGRATOR_
#define _FETIMEINTEGRATOR_

#include "feNumber.h"
#include "feSolution.h"
#include "feMesh.h"

class feTimeIntegrator{

protected:
  int _nDofs; // NDDL
  int _nSol;
  std::vector<double> _t; // TEMPS
  std::vector<double> _sol; // U
  std::vector<double> _fResidual; // F
  std::vector<double> _cn;
  std::vector<double> _d;

  // Donner aussi la solution ?

public:
	feTimeIntegrator(int nSol, double tn, feMetaNumber *metaNumber)
    : _nSol(nSol), _nDofs(metaNumber->getNbDOFs())
	{
    _t.resize(_nSol);
    _t[0] = tn;
    _sol.resize(_nSol*_nDofs);
    _fResidual.resize(_nSol*_nDofs);
    _cn.resize(3);
    _d.resize(_nDofs);
	};
	~feTimeIntegrator() {
  }

  void initialize(feSolution *sol, feMesh *mesh, feMetaNumber *metaNumber);

};

#endif