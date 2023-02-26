#ifndef _FENEWRECOVERY_
#define _FENEWRECOVERY_

#include "feMesh.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feVertex.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feFunction.h"
#include "feMatrixInterface.h"

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

#include <iostream>
#include <fstream>

#if defined(HAVE_EIGEN)
#include "../contrib/Eigen/Dense"
#endif

enum class PPR {RECOVERY, DERIVATIVE};

// Patches of elements around vertices of a single connectivity
class feNewPatch
{
protected:
  // Number of connectivity vertices = number of patches
  int _nVertices;
  // Global indices in mesh->_vertices
  std::vector<int> _vertices;
  // Number of layers of elements around the vertex
  int _nLayers;

  int _nNodePerElm; // !!! Modifier pour nVert pour des P2+ geometriques ?
  int _nEdgePerElm;

  std::map<int, std::set<int> > vertToElems; // A set of element tags for each mesh vertex
  std::map<int, std::set<int> > edgeToElems; // A set of element tags for each mesh edge

public:
  feNewPatch(const feCncGeo *cnc, feMesh *mesh);
  ~feNewPatch() {}

  std::vector<int> &getVertices() { return _vertices; }
  std::set<int> &getPatch(int iVertex) { return vertToElems[iVertex]; }
  std::set<int> &getEdgePatch(int iEdge) { return edgeToElems[iEdge]; }
};

typedef std::map<int, Eigen::MatrixXd > matMap;
typedef std::map<int, SquareMatrix > matMap2;

class feNewRecovery
{
public:
  feCncGeo *_cnc;

  int _nElm;
  int _nNodePerElm;
  std::vector<feInt> _adr;
  std::vector<double> _solution;

  feMesh *_mesh;
  feSolution *_sol;

  feSpace *_intSpace;
  feSpace *_geoSpace;

  feNewPatch *_patch;

  std::map<int, std::map<int, std::vector<double> > > recoveryCoeff; // #vert : {#rec , coeffs}

  std::map<int, std::map<int, std::map<int, std::vector<double> > > >
    recoveryCoeffOnEdges; // #edge : {#dof : {#rec , coeffs}}

  std::map<int, std::map<int, std::vector<double> > > derivativeCoeff; // #vert : {#der , coeffs}

  std::map<int, std::map<int, std::map<int, std::vector<double> > > >
    derivativeCoeffOnEdges; // #edge : {#dof : {#rec , coeffs}}

  std::map<int, std::vector<double> > errorCoeff; // #vert : coeffs

  /* All recovered functions (solution and derivatives) evaluated at the vertices dof : #recovery :
    {#vertex : val} They are pushed after they are created : 1D : u, dudx, d2udx2, d3udx3, ...
      2D : u, dudx, dudy, d2udx2, d2udxy, d2udyx, d2udy2, d3udx3, ...  */
  std::vector<std::vector<double> > derivAtVertices;
  // All recovered functions evaluated at the edges dof : #derivative : {#edge-1 : val}
  std::vector<std::vector<double> > derivAtEdges;
protected:
  int _degSol; // The degree of the FE solution : degree of recovery = deg of the FE solution +1
  int _degRec;

  int _dim;
  int _nTotalRecoveries;
  int _nTotalDerivations;
  int _dimRecovery;
  int _dimDerivation;
  int _dim2Derivation;

  std::vector<int> _expX_recov;
  std::vector<int> _expY_recov;
  std::vector<int> _expX_deriv;
  std::vector<int> _expY_deriv;

  // To interpolate the recovery/derivatives
  std::vector<double> _FIELD;
  std::vector<double> _geoCoord;
  std::vector<double> _pos;
  std::vector<double> _basisRecovery;
  std::vector<double> _basisDerivative;

  matMap _inverseMassMatrices;
  matMap2 _inverseMassSquareMatrices;

public:
  feNewRecovery(feSpace *space, feMesh *mesh, feSolution *sol,
             std::string meshName = "",
             std::string metricMeshName = "",
             bool append = false);
  feNewRecovery(feSpace *space, feMesh *mesh, std::string recoveryFile);
  ~feNewRecovery() { delete _patch; }

  int getDim() { return _dim; }
  int getDimRecovery() { return _dimRecovery; }
  int getDegreeSolution() { return _degSol; }

  std::vector<int> &getVertices() { return _patch->getVertices(); }

  // std::vector<int> &getXExponentsRecovery() { return _expXRecovery; }
  // std::vector<int> &getYExponentsRecovery() { return _expYRecovery; }
  // std::vector<int> &getZExponentsRecovery() { return _expZRecovery; }
  // std::map<int, std::vector<double> > &getErrorCoefficients() { return errorCoeff; }

  // void writeRecovery(std::string fileName);

  double evaluateRecovery(PPR recoveredField,
    const int index, const double *x);
  double evaluateRecoveryAtVertex(PPR recoveredField,
    const int index, const int vertex);
  double evaluateRecoveryAtQuadNode(PPR recoveredField,
    const int index, const int iElm, const int iQuadNode);

private:
  void setDimensions();
  void setPolynomialExponents();
  void computeVertexMassMatrices2D();
  void inverseMassMatrices();
  void computeRHSAndSolve(int numRecoveries, int nRecoveredFields, int nRecoveredDerivatives, int iDerivative);
  void computeDerivative(int nRecoveredFields, int indexRecovery, int iDerivative, std::ostream &output);
  double evaluatePolynomial(PPR recoveredField,
    const int index, const int vertex, const double *xLoc);
  void computeHomogeneousErrorPolynomials();
  void computeRecoveryAtAllElementDOF(PPR recoveredField, const int index, const int iElm);
};

#endif
