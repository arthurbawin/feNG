#ifndef _FECNCGEO_
#define _FECNCGEO_

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>

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
  // Sequential unique number
  int _tag;
  // Dimension
  int _dim;
  // Elements geometry (see enum above)
  geometryType _geometry;
  // Interpolation functions used to describe the geometry (see enum above)
  geometricInterpolant _interpolant;

  // Number of vertices in the Physical Entity
  int _nVertices;
  // Number of vertices per element
  int _nVerticesPerElm;
  // Number of elements
  int _nElements;
  // Number of edges
  int _nEdges;
  // Number of facets (edges for now) per element
  int _nEdgesPerElem;

  // Global numbering of the entities (vertex/edge/element)
  // on this connectivity.

  // The global tag of vertices classified on this connectivity: size = nVertices
  std::vector<int> _connecVerticesOnly;
  // The global tag of vertices of elements on this connectivity: size = nElm * nVerticesPerElm
  std::vector<int> _connecVertices;
  // The global tag of edges classified on this connectivity: size = nEdges
  std::vector<int> _connecEdgesOnly;
  // The global tag of edges of elements on this connectivity: size = nElm * nEdgesPerElm
  std::vector<int> _connecEdges;
  std::vector<int> _connecFaces;
  std::vector<int> _connecElem;

  // Pointers to the FE space used to interpolate the geometry
  // and to the mesh
  feSpace *_geometricInterpolant;
  feMesh *_mesh;

  // Determinants of the jacobian matrix of the reference-to-physical
  // transformation, stored at all quadrature nodes of all elements
  // on this connectivity. Size = nQuad x nElm
  std::vector<double> _J;

  // Mesh coloring
  std::vector<int> _nbElmPerNode;
  std::vector<int> _nbElmPerElm;
  std::vector<std::vector<int> > _listElmPerNode;
  std::vector<std::vector<int> > _listElmPerElm;

  int _nbColor;
  std::vector<int> _elmToColor;
  std::vector<int> _nbElmPerColor;
  std::vector<std::vector<int> > _listElmPerColor;

public:
  struct coloring {
    int numColors;
    std::vector<int> elem2Color;
    std::vector<int> numElemPerColor;
    std::vector<std::vector<int> > elementsInColor;
  };

  coloring _coloring;

public:
  // Create a geometric connectivity. Called when parsing the mesh.
  feCncGeo(const int tag, const int dimension, const int nVerticesPerElement, const int nElements,
           const int nEdgesPerElement, const std::string &ID, const geometryType geometry,
           const geometricInterpolant interpolant, feSpace *space, std::vector<int> connecVertices,
           std::vector<int> connecElem = std::vector<int>(),
           std::vector<int> connecEdges = std::vector<int>(),
           std::vector<int> connecFaces = std::vector<int>());
  ~feCncGeo() {}

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

  // Set quadrature rule for the interpolation space.
  // Triggers the computation of jacobian determinants.
  feStatus setQuadratureRule(feQuadrature *rule);

  void setMeshPtr(feMesh *mesh) { _mesh = mesh; }

  feStatus computeJacobians();
  const std::vector<double> &getJacobians() const { return _J; }

  void computeElementTransformation(std::vector<double> &elementCoord, const int iQuadNode,
                                    const double jac, ElementTransformation &transformation) const;

  const std::vector<int> &getVerticesConnectivity() const { return _connecVertices; }
  const std::vector<int> &getEdgeConnectivity() const { return _connecEdges; }
  const std::vector<int> &getElemConnectivity() const { return _connecElem; }

  // Get and set vertex connectivity
  int getUniqueVertexConnectivity(const int iVertex) const;
  int getVertexConnectivity(const int iVertex) const;
  int getVertexConnectivity(const int numElem, const int iVertex) const;
  void setVertexConnectivity(const int numElem, const int iVertex, const int val);

  // Get and set element connectivity
  int getElementConnectivity(const int numElem) const;
  void setElementConnectivity(const int numElem, const int val);

  // Get and set edge connectivity
  int getUniqueEdgeConnectivity(const int iEdge) const;
  int getEdgeConnectivity(const int numElem, const int iEdge) const;
  void setEdgeConnectivity(const int numElem, const int iEdge, const int val);

  void colorElements(int coloringAlgorithm);

  int getNbColor() const { return _nbColor; };
  const std::vector<int> &getColorElm() const { return _elmToColor; };
  const std::vector<int> &getNbElmPerColor() const { return _nbElmPerColor; };
  const std::vector<std::vector<int> > &getListElmPerColor() const { return _listElmPerColor; };
  int getNbElmPerColorI(int i) const { return _nbElmPerColor[i]; };
  const std::vector<int> &getListElmPerColorI(int i) const { return _listElmPerColor[i]; };
  int getElmColored(int iColor, int iElmC) const { return _listElmPerColor[iColor][iElmC]; };
};

#endif
