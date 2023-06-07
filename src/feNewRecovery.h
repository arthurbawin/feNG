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

enum class PPR { RECOVERY, DERIVATIVE };

class feMetric;

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
  std::map<int, std::set<int> > vertToVerts; // A set of vertices tags for mesh vertices (the vertices of its elem patch)
  std::map<int, std::set<int> > vertToInteriorVerts; // A set of interior vertices tags for boundary mesh vertices (the interior vertices of its elem patch)
  std::map<int, std::set<int> > edgeToElems; // A set of element tags for each mesh edge
  std::unordered_map<int, bool> isOnEdge;

  std::map<int, bool> _isVertexBoundary;
  std::map<int, std::pair<double, double>> _scaling;

public:
  feNewPatch(const feCncGeo *cnc, feMesh *mesh, bool reconstructAtHighOrderNodes, feMetric *metricField = nullptr);
  ~feNewPatch() {}

  std::vector<int> &getVertices() { return _vertices; }
  std::map<int, std::set<int>> &getElementPatches() { return vertToElems; };
  std::map<int, std::set<int>> &getVerticesPatches() { return vertToVerts; };
  std::set<int> &getPatch(int iVertex) { return vertToElems[iVertex]; }
  std::set<int> &getPatchOfInteriorVertices(int iVertex) { return vertToInteriorVerts.at(iVertex); }
  std::set<int> &getEdgePatch(int iEdge) { return edgeToElems[iEdge]; }
  std::map<int, bool> &getBoundaryVertices() { return _isVertexBoundary; };
  std::map<int, std::pair<double, double>> &getScaling() { return _scaling; };

  bool isVertexBoundary(int iVertex) { return _isVertexBoundary.at(iVertex); };
  bool isVertexOnEdge(int iVertex) { return isOnEdge.at(iVertex); };
};

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
  feMetaNumber *_numbering;

  feSpace *_intSpace;
  feSpace *_geoSpace;

  feNewPatch *_patch;

  std::map<int, std::map<int, std::vector<double> > > recoveryCoeff; // #vert : {#rec , coeffs}

  // All recoveries' independant term condensed in a vector
  // The ith recovery at vertex vertex is at position recoveryIndependantTerm[ numTotalRecoveries * vertex + i ]
  std::vector<double> recoveryIndependantTerm; 

  std::map<int, std::map<int, std::map<int, std::vector<double> > > >
    recoveryCoeffOnEdges; // #edge : {#dof : {#rec , coeffs}}

  std::map<int, std::map<int, std::vector<double> > > derivativeCoeff; // #vert : {#der , coeffs}

  // All derivatives' independant term condensed in a vector
  // The ith recovery at vertex vertex is at position derivativeIndependantTerm[ numTotalRecoveries * vertex + i ]
  std::vector<double> derivativeIndependantTerm;

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

  int _numTotalRecoveries;
  int _numTotalDerivatives;

  int _dimRecovery;
  int _dimDerivation;
  int _dim2Derivation;

  bool _reconstructAtHighOrderNodes;

  std::vector<int> _expX_recov;
  std::vector<int> _expY_recov;
  std::vector<int> _expX_deriv;
  std::vector<int> _expY_deriv;

  // 1D quadrature rule to compute transformed edges lengths
  int _nQuad1d;
  std::vector<double> _wQuad1d;
  std::vector<double> _xQuad1d;

  // To interpolate the recovery/derivatives
  // std::vector<double> _FIELD;
  // std::vector<double> _geoCoord;
  // std::vector<double> _pos;

  std::map<int, Eigen::MatrixXd> _inverseMassMatrices;
  std::map<int, SquareMatrix> _inverseMassSquareMatrices;
  std::map<int, Matrix> _leastSquaresMatrices;

public:
  feNewRecovery(feSpace *space,
                feMesh *mesh,
                feSolution *sol,
                std::string meshName = "",
                std::string metricMeshName = "",
                bool reconstructAtHighOrderNodes = false,
                bool append = false,
                feMetric *metricField = nullptr,
                feMetaNumber *numbering = nullptr);
  feNewRecovery(feSpace *space, feMesh *mesh, std::string recoveryFile);
  ~feNewRecovery() { delete _patch; }

  int getDim() { return _dim; }
  int getDimRecovery() { return _dimRecovery; }
  int getDegreeSolution() { return _degSol; }

  std::vector<int> &getVertices() { return _patch->getVertices(); }

  // std::vector<int> &getXExponentsRecovery() { return _expXRecovery; }
  // std::vector<int> &getYExponentsRecovery() { return _expYRecovery; }
  // std::vector<int> &getZExponentsRecovery() { return _expZRecovery; }
  std::map<int, std::vector<double> > &getErrorCoefficients() { return errorCoeff; }

  // void writeRecovery(std::string fileName);

  double evaluateRecovery(PPR recoveredField, const int index, const double *x, bool averageEvaluations = false);
  double evaluateRecovery(PPR recoveredField, const int index, const double *x, const std::vector<int> &elementsToSearch);
  double evaluateRecoveryLinear(PPR recoveredField, const int index, const double *x);
  double evaluateRecoveryAtVertex(PPR recoveredField, const int index, const int vertex);
  double evaluateRecoveryAtQuadNode(PPR recoveredField, const int index, const int iElm,
                                    const int iQuadNode);
  double evaluateRecoveryAtQuadNodeLinear(PPR recoveredField, const int index, const int iElm,
                                    const int iQuadNode);

  // Recompute recovery and gradient during curving
  void recomputeLeastSquareMatrix(const int vertex);
  void recomputeRHSAndSolve(const int vertex);
  void recomputeFirstDerivative(const int vertex);

private:
  void setDimensions();
  void setPolynomialExponents();
  void computeVertexMassMatrices2D();
  void computeVertexLeastSquareMatrices2D();
  void inverseMassMatrices();
  void computeRHSAndSolve(int numRecoveries, int nRecoveredFields, int nRecoveredDerivatives,
                          int iDerivative);
  void computeRHSAndSolve_noIntegral(int numRecoveries, int nRecoveredFields, int nRecoveredDerivatives,
                          int iDerivative);
  void computeDerivative(int nRecoveredFields, int indexRecovery, int iDerivative,
                         std::ostream &output);
  double evaluatePolynomial(PPR recoveredField, const int index, const int vertex,
                            const double *xLoc);
  double evaluatePolynomialAtBoundaryVertex(PPR recoveredField, const int index, const int vertex,
                                            const std::set<int> &connectedInteriorVertices);
  double averageRecoveryEvaluations(PPR recoveredField, const int index, const int iElm, double *r);
  double averageRecoveryEvaluationsAtQuadNode(PPR recoveredField, const int index,
                                                 const int iElm, const int iQuadNode);
  void computeHomogeneousErrorPolynomials();
  void computeTransformedEdgeLength(std::vector<double> &geoCoord,
                                    std::vector<double> &dxdr,
                                    std::vector<double> &dxds,
                                    int whichEdge, double &firstHalf, double &secondHalf);
  void computeRecoveryAtAllElementDOF(PPR recoveredField,
                                      const int index,
                                      const int iElm,
                                      std::vector<double> &recoveryAtDOFS,
                                      bool computeAtP1VerticesOnly = false);
  void computeRecoveryAtAllElementDOF2(PPR recoveredField,
                                      const int index,
                                      const int iElm,
                                      std::vector<double> &recoveryAtDOFS);
};

#endif
