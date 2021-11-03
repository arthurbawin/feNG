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

#if defined(HAVE_EIGEN)
#include "../contrib/Eigen/Dense"
#endif

// A structure for the patches of elements around vertices of a connectivity
class fePatch {
protected:
  int _nVertices; // Number of mesh vertices = number of patches
  std::vector<int> _vertices; // The indices in mesh->_vertices

  int _nNodePerElm; // !!! Modifier pour nVert pour des P2+ geometriques ?
  int _nEdgePerElm;

  std::map<int, std::set<int> > vertToElems; // A set of element tags for each mesh vertex
  std::map<int, std::set<int> > edgeToElems; // A set of element tags for each mesh edge

public:
  fePatch(feCncGeo *cnc, feMesh *mesh);
  ~fePatch() {}

  std::vector<int> &getVertices() { return _vertices; }
  std::set<int> &getPatch(int iVertex) { return vertToElems[iVertex]; }
  std::set<int> &getEdgePatch(int iEdge) { return edgeToElems[iEdge]; }
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

#if defined(HAVE_EIGEN)
  std::map<int, Eigen::Matrix<double, 6, 6> > lsInvAtVertices6;
  std::map<int, Eigen::Matrix<double, 6, 6> > lsInvAtVertices6OnEdges;
  std::map<int, Eigen::Matrix<double, 10, 10> > lsInvAtVertices10;
  std::map<int, Eigen::Matrix<double, 10, 10> > lsInvAtVertices10OnEdges;
#endif

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

  std::map<int, std::map<int, std::vector<double> > > recoveryCoeff; // #vert : {#rec , coeffs}
  std::map<int, std::map<int, std::map<int, std::vector<double> > > >
    recoveryCoeffOnEdges; // #edge : {#dof : {#rec , coeffs}}
  std::map<int, std::map<int, std::vector<double> > > derivativeCoeff; // #vert : {#der , coeffs}
  std::map<int, std::map<int, std::map<int, std::vector<double> > > >
    derivativeCoeffOnEdges; // #edge : {#dof : {#rec , coeffs}}
  std::map<int, std::vector<double> > errorCoeff; // #vert : coeffs

public:
  feRecovery(feMetaNumber *metaNumber, feSpace *space, feMesh *mesh, feSolution *sol,
             std::vector<double> &norm, feFunction *solRef, std::string meshName = "",
             std::string metricMeshName = "", feVectorFunction *solRefGrad = nullptr,
             feVectorFunction *solRefHess = nullptr, bool append = false);
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
  void matrixInverseEigen();
  // void freeStructures();
  // void solveLeastSquare(int indRecovery, bool recoverDerivative);
  void solveLeastSquareEigen(int indRecovery, bool recoverDerivative);
  // void derivative(int indRecovery, int iDerivative, FILE* derivativeFile);
  void derivative(int indRecovery, int iDerivative, std::ostream &output);
  void getErrorPolynomials();
  void estimateError(std::vector<double> &norm, feFunction *solRef);
  void estimateH1Error(std::vector<double> &norm, feVectorFunction *solRefGrad);
  void estimateDudxError(std::vector<double> &norm, feVectorFunction *solRefGrad);
  void estimateHessError(std::vector<double> &norm, feVectorFunction *solRefHess);
};

#endif