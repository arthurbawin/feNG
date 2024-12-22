#ifndef _FECNCGEO_
#define _FECNCGEO_

#include <memory>
#include "feColoring.h"
#include "feMessage.h"

// Supported geometries are only POINT, LINE and TRI so far
enum class geometryType { POINT, LINE, TRI, QUAD, TET, HEX };

// Supported geometric interpolants (Lagrange interpolation only)
enum class geometricInterpolant {
  NONE,
  POINTP0,
  LINEP1,
  LINEP2,
  LINEP3,
  TRIP1,
  TRIP2,
  TRIP3,
  QUADP1,
  QUADP2,
  TETP1,
  TETP2
};

std::string toString(geometryType t);
std::string toString(geometricInterpolant t);
int getGeometricInterpolantDegree(geometricInterpolant t);

// Reference-to-physical transformation jacobian matrix dx(r)/dr
// and its inverse dr(x)/dx.
typedef struct TransformationStruct {
  double jac;
  double dxdr[3];
  double dxds[3];
  double dxdt[3];
  double drdx[3];
  double drdy[3];
  double drdz[3];
} ElementTransformation;

class feSpace;
class feMesh;
class feQuadrature;

//
// A geometric connectivity is associated to each Physical Entities
// in the mesh.
//
class feCncGeo
{
protected:
  // Name of the Physical Entity
  std::string _ID;
  // Unique sequential number
  int _tag;
  // Dimension of this Physical entity and dimension of the ambient space
  int _dim;
  int _ambientDim;
  // Elements geometry
  geometryType _geometry;
  // Interpolation functions used to describe the geometry (see enum above)
  geometricInterpolant _interpolant;

  // Pointers to the mesh and to the FE space used to interpolate the geometry
  feMesh *_mesh;
  feSpace *_geometricInterpolant;

  // Number of elements of the dimension of the Physical Entity
  // (i.e., elements of highest dimension for this entity)
  int _nElements;

  // Total number of mesh elements classified on this Physical Entity
  int _nVertices;
  int _nEdges;
  int _nTriangles;
  int _nTetrahedra;

  // Number of lower-dimensional elements for the elements of this entity
  int _nVerticesPerElm;
  int _nEdgesPerElem;
  int _nTrianglesPerElem;

  //
  // Global numbering of the mesh elements on this connectivity.
  //
  std::vector<int> _connecElem;

  // The global tags of mesh elements classified on this connectivity
  // Size : _nVertices // _nEdges // _nTriangles // _nTetrahedra
  std::vector<int> _connecVerticesOnly;
  std::vector<int> _connecEdgesOnly;
  std::vector<int> _connecTriOnly;
  std::vector<int> _connecTetOnly;

  // The global tags of boundary elements on this connectivity
  // Size : _nElements * (_nVerticesPerElem // _nEdgesPerElem // _nTrianglesPerElem)
  std::vector<int> _connecVertices;
  std::vector<int> _connecEdges;
  std::vector<int> _connecTriangles;

  // Computing the Jacobian determinants on a Physical Entity
  // requires knowing the quadrature points at which they are evaluated,
  // thus setting a quadrature rule, which requires knowing the maximum
  // polynomial degree to integrate exactly, which is not required to create a mesh.
  // The quadrature rule is set whenever a finite element
  // space is created on this entity.
  // Alternatively, we could:
  // - require a polynomial degree to create the mesh structure (a bit counterintuitive)
  // - create a rule with an arbitrary degree
  // - ... ?
  bool _jacobiansWereComputed = false;

  // Determinants of the jacobian matrix of the reference-to-physical
  // transformation, stored at all quadrature nodes of all elements
  // on this connectivity. Size = nQuad x nElm
  std::vector<double> _J;
  std::vector<double> _elementsVolume;
  std::vector<double> _minimumScaledJacobianControlCoeffs;

  // ======= To be replaced by just a feColoring ==============
  // Mesh coloring
  std::vector<int> _nbElmPerNode;
  std::vector<int> _nbElmPerElm;
  std::vector<std::vector<int> > _listElmPerNode;
  std::vector<std::vector<int> > _listElmPerElm;

  int _nbColor;
  std::vector<int> _elmToColor;
  std::vector<int> _nbElmPerColor;
  std::vector<std::vector<int> > _listElmPerColor;
  // ==========================================================

public:
  std::unique_ptr<feColoring> _mycoloring;

public:
  // Create a geometric connectivity. Called when parsing the mesh.
  feCncGeo(const int tag, const int dimension, const int ambientDimension,
           const int nVerticesPerElement,
           const int nElements,
           const int nEdgesPerElement, 
           const int nTrianglesPerElement,
           const std::string &ID, const geometryType geometry,
           const geometricInterpolant interpolant, feSpace *space,
           std::vector<int> connecVertices,
           std::vector<int> connecElem = std::vector<int>(),
           std::vector<int> connecEdges = std::vector<int>(),
           std::vector<int> connecTriangles = std::vector<int>());
  ~feCncGeo() {
    // delete _mycoloring;
  }

  const std::string &getID() const { return _ID; }
  int getTag() const { return _tag; }
  int getDim() const { return _dim; }
  int getNumVertices() const { return _nVertices; }
  int getNumVerticesPerElem() const { return _nVerticesPerElm; }
  int getNumElements() const { return _nElements; }
  int getNumEdges() const { return _nEdges; }
  int getNumEdgesPerElem() const { return _nEdgesPerElem; }
  geometryType getGeometry() const { return _geometry; };
  geometricInterpolant getInterpolant() const { return _interpolant; };
  feSpace *getFeSpace() const { return _geometricInterpolant; }
  feMesh *getMeshPtr() { return _mesh; }
  void setMeshPtr(feMesh *mesh) { _mesh = mesh; }

  const std::vector<int> &getVerticesConnectivity() const { return _connecVertices; }
  const std::vector<int> &getEdgeConnectivity() const { return _connecEdges; }
  const std::vector<int> &getElemConnectivity() const { return _connecElem; }

  //
  // Get (and rarely set) connectivities on this Physical Entity
  //
  int getUniqueVertexConnectivity(const int iVertex) const;
  int getVertexConnectivity(const int iVertex) const;
  int getVertexConnectivity(const int numElem, const int iVertex) const;
  // void setVertexConnectivity(const int numElem, const int iVertex, const int val);
  int getElementConnectivity(const int numElem) const;
  void setElementConnectivity(const int numElem, const int val);
  int getUniqueEdgeConnectivity(const int iEdge) const;
  int getEdgeConnectivity(const int numElem, const int iEdge) const;
  // void setEdgeConnectivity(const int numElem, const int iEdge, const int val);
  int getFaceConnectivity(const int numElem, const int iFace) const;

  const std::vector<double> &getJacobians() const { return _J; }
  const std::vector<double> &getElementsVolume() const
    { return _elementsVolume; };
  const std::vector<double> &getMinimumScaledJacobianControlCoeffs() const
    { return _minimumScaledJacobianControlCoeffs; };

  bool jacobiansWereComputed() const { return _jacobiansWereComputed; };
  feStatus computeJacobians(const bool ignoreNegativeJacobianWarning = false);
  feStatus recomputeElementJacobian(const int iElm);

  // Set quadrature rule for the interpolation space.
  // Triggers the computation of jacobian determinants.
  feStatus setQuadratureRule(feQuadrature *rule);
  feStatus computeMinimumScaledJacobianControlCoefficients();

  void computeElementTransformation(std::vector<double> &elementCoord, const int iQuadNode,
                                    const double jac, ElementTransformation &transformation) const;

  feStatus computeNormalVectors(std::vector<double> &normalVectors) const;

  // ========== To wrap in a feColoring ================
  void colorElements(int coloringAlgorithm);
  int getNbColor() const { return _nbColor; };
  const std::vector<int> &getColorElm() const { return _elmToColor; };
  const std::vector<int> &getNbElmPerColor() const { return _nbElmPerColor; };
  const std::vector<std::vector<int> > &getListElmPerColor() const { return _listElmPerColor; };
  int getNbElmPerColorI(int i) const { return _nbElmPerColor[i]; };
  const std::vector<int> &getListElmPerColorI(int i) const { return _listElmPerColor[i]; };
  int getElmColored(int iColor, int iElmC) const { return _listElmPerColor[iColor][iElmC]; };
  void printColoringStatistics();
  void printColoring(std::string fileName);
  void printColoring2(std::string fileName);
  // ===================================================

  // Write constant field on an element to a .pos file (Gmsh visualization)
  void writeElementToPOS(FILE *posFile, const std::vector<double> &elementCoord, const double value) const;
};

#endif
