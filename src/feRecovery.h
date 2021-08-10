#ifndef _FERECOVERY_
#define _FERECOVERY_

#include "feMesh.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feSpaceTriangle.h"
#include "feVertex.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feFunction.h"

#include "petscksp.h"

#include <iostream>
#include <fstream>

// A structure for the patches of elements around vertices of a connectivity
class fePatch {
protected:
  int _nVertices; // Number of mesh vertices = number of patches
  std::vector<int> _vertices; // The indices in mesh->_vertices

  int _nNodePerElm; // !!! Modifier pour nVert pour des P2+ geometriques ?

  std::map<int, std::set<int> > vertToElems; // A set of element tags for each mesh vertex

public:
  fePatch(feCncGeo *cnc);
  ~fePatch() {}

  std::vector<int> &getVertices() { return _vertices; }
  std::set<int> &getPatch(int iVertex) { return vertToElems[iVertex]; }
};

class feRecovery {
protected:
  int _degSol; // The degree of the FE solution : degree of recovery = deg of the FE solution +1

  int _dim;
  int _nTotalRecoveries;
  int _nTotalDerivations;
  int _dimRecovery;
  int _dimDerivation;

  std::vector<int> _expXRecovery;
  std::vector<int> _expYRecovery;
  std::vector<int> _expZRecovery;
  std::vector<int> _expX;
  std::vector<int> _expY;
  std::vector<int> _expZ;

  Mat A;
  Vec b, c, res;
  KSP ksp;
  PC preconditioner;

  PetscInt *indI;
  PetscInt *indJ;
  PetscScalar *valA;
  PetscScalar *valb;

  std::map<int, std::map<int, std::vector<double> > > recoveryCoeff; // #vert : {#rec , coeffs}
  std::map<int, std::map<int, std::vector<double> > > derivativeCoeff; // #vert : {#der , coeffs}
  std::map<int, std::vector<double> > errorCoeff; // #vert : coeffs

public:
  int _nElm;
  int _nNodePerElm;

  feMetaNumber *_metaNumber;
  feMesh *_mesh;
  feSolution *_sol;

  feSpace *_intSpace;
  feSpace *_geoSpace;

  fePatch *_patch;
  feCncGeo *_cnc;

public:
  feRecovery(feMetaNumber *metaNumber, feSpace *space, feMesh *mesh, feSolution *sol,
             std::vector<double> &norm, feFunction *solRef, std::string meshName = "", std::string metricMeshName = "",
             feVectorFunction *solRefGrad = nullptr);
  ~feRecovery() { delete _patch; }

  int getDim() { return _dim; }
  int getDimRecovery() { return _dimRecovery; }
  int getDegreeSolution() { return _degSol; }

  std::vector<int> &getVertices() { return _patch->getVertices(); }

  std::vector<int> &getXExponentsRecovery() { return _expXRecovery; }
  std::vector<int> &getYExponentsRecovery() { return _expYRecovery; }
  std::vector<int> &getZExponentsRecovery() { return _expZRecovery; }
  std::map<int, std::vector<double> > &getErrorCoefficients() { return errorCoeff; }

  void allocateStructures();
  void freeStructures();
  void solveLeastSquare(int indRecovery, bool recoverDerivative);
  // void derivative(int indRecovery, int iDerivative, FILE* derivativeFile);
  void derivative(int indRecovery, int iDerivative, std::ostream &output);
  void getErrorPolynomials();
  void estimateError(std::vector<double> &norm, feFunction *solRef);
  void estimateH1Error(std::vector<double> &norm, feVectorFunction *solRefGrad);
};

#endif