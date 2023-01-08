#ifndef _FEBILINEARFORM_
#define _FEBILINEARFORM_

#include "feNG.h"
#include "feMesh.h"
#include "feSpace.h"
#include "feSysElm.h"
#include "feQuadrature.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feCncGeo.h"

class feBilinearForm
{
protected:
  feSysElm *_sysElm;

public:
  // These members are used to compute the elementary system.
  // Since friendship is not inherited by feSysElm derived class,
  // they are public for now although it's probably not optimal.
  std::vector<feSpace *> _intSpace;
  feCncGeo *_cnc;
  feSpace *_geoSpace;
  std::vector<double> _geoCoord;

  double _c0; // First coefficient of the BDF expansion
  double _tn; // Current time
  double _dt; // Time step

  int _numElem;

  double **_Ae;
  double *_Be;

protected:
  int _cncGeoTag;
  std::string _cncGeoID;
  int _nCoord;
  int _nGeoNodes;
  int _nGeoElm;

  int _nQuad;
  int _degQuad;
  std::vector<double> _w;
  std::vector<double> _x;
  std::vector<double> y;
  std::vector<double> z;

  std::vector<int> _iVar;
  std::vector<int> _jVar;
  feInt _niElm;
  feInt _njElm;
  std::vector<feInt> _adrI;
  std::vector<feInt> _adrJ;  

public:
  std::vector<std::vector<feInt>> _adr;

  // The solution at DOFs on the current element for all FE spaces
  std::vector<std::vector<double>> _sol;
  std::vector<std::vector<double>> _solDot;

  // The solution on the previous and next element (for e.g. DG fluxes)
  std::vector<std::vector<double>> _solPrev;
  std::vector<std::vector<double>> _solNext;

public:
  // To compute the elementary matrix with finite differences
  double *_R0;
  double *_Rh;
  double _h0;

  void (feBilinearForm::*ptrComputeMatrix)(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                           int numElem);

public:
  feBilinearForm(std::vector<feSpace *> space, feMesh *mesh, int degQuad, feSysElm *sysElm);
  feBilinearForm(const feBilinearForm &f);
  ~feBilinearForm();

  feCncGeo *getCncGeo() { return _cnc; }
  int getCncGeoTag() { return _cncGeoTag; }

  bool hasMatrix() { return _sysElm->hasMatrix(); }

  feInt getNiElm() { return _niElm; }
  feInt getNjElm() { return _njElm; }
  std::vector<feInt> &getAdrI() { return _adrI; }
  std::vector<feInt> &getAdrJ() { return _adrJ; }

  double **getAe() { return _Ae; }
  double *getBe() { return _Be; }

  elementSystemType getID() { return _sysElm->getID(); }
  std::string getIDName() { return _sysElm->getIDName(); }

  void initialize_vadij_only(feMetaNumber *metaNumber, int numElem);
  void initialize(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);

  void computeMatrix(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);
  void computeMatrixAnalytical(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                               int numElem);
  void computeMatrixFiniteDifference(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol,
                                     int numElem);
  void computeResidual(feMetaNumber *metaNumber, feMesh *mesh, feSolution *sol, int numElem);

  double getMatrixNorm();
  double getResidualNorm();
};

#endif