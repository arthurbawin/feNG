#ifndef _FECNCGEO_
#define _FECNCGEO_

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>

#include "feMessage.h"

// Supported geometries are only POINT, LINE and TRI so far
// typedef enum { POINT, LINE, TRI, QUAD, TET, HEX } geometryType;

enum class geometryType {
  POINT, LINE, TRI, QUAD, TET, HEX
};

// Supported geometric interpolants (Lagrange interpolation only)
// Somewhat redundant
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

inline std::string toString(geometryType t)
{
  switch(t) {
    case geometryType::POINT: return "POINT";
    case geometryType::LINE: return "LINE";
    case geometryType::TRI: return "TRI";
    case geometryType::QUAD: return "QUAD";
    case geometryType::TET: return "TET";
    case geometryType::HEX: return "HEX";
    default: return "Unknown geometry";
  }
}

inline std::string toString(geometricInterpolant t)
{
  switch(t) {
    case geometricInterpolant::NONE:
      return "NONE";
    case geometricInterpolant::POINTP0:
      return "POINTP0";
    case geometricInterpolant::LINEP1:
      return "LINEP1";
    case geometricInterpolant::LINEP2:
      return "LINEP2";
    case geometricInterpolant::LINEP3:
      return "LINEP3";
    case geometricInterpolant::TRIP1:
      return "TRIP1";
    case geometricInterpolant::TRIP2:
      return "TRIP2";
    case geometricInterpolant::TRIP3:
      return "TRIP3";
    default:
      return "Unknown geometric interpolant";
  }
}

class feSpace;
class feMesh;
class feQuadrature;

class feCncGeo
{
protected:
  std::string _ID;
  int _tag;
  int _dim;

  geometryType _geometry;

  geometricInterpolant _interpolant;

  int _nNod;
  int _nNodPerElm; // Nombre de noeuds par element
  int _nElm; // Nombre d'elements dans cette connectivite
  int _nEdg;
  std::vector<int> _connecNodes; // Connectivite des noeuds
  std::vector<int> _connecElem; // Connectivite des elements
  std::vector<int> _connecEdges; // Connectivite des aretes
  std::vector<int> _connecFaces; // Connectivite des faces

  feSpace *_space;
  feMesh *_mesh;

  std::vector<double> _J; // Jacobiens

  // For elements coloring
  std::vector<int> _nbElmPerNode;
  std::vector<int> _nbElmPerElm;
  std::vector<std::vector<int> > _listElmPerNode;
  std::vector<std::vector<int> > _listElmPerElm;

  int _nbColor;
  std::vector<int> _elmToColor;
  std::vector<int> _nbElmPerColor;
  std::vector<std::vector<int> > _listElmPerColor;

public:
  feCncGeo(int tag, int dim, int nNod, int nElm, int nEdg, std::string ID, 
    geometryType geometry, geometricInterpolant interpolant,
           feSpace *space, std::vector<int> connecNodes,
           std::vector<int> connecElem = std::vector<int>(),
           std::vector<int> connecEdges = std::vector<int>(),
           std::vector<int> connecFaces = std::vector<int>())
    : _ID(ID), _tag(tag), _dim(dim), _geometry(geometry), _interpolant(interpolant), _nNodPerElm(nNod), _nElm(nElm), _nEdg(nEdg),
      _connecNodes(connecNodes), _connecElem(connecElem), _connecEdges(connecEdges),
      _connecFaces(connecFaces), _space(space)
  {
    if(connecElem.size() == 0)
      _connecElem.resize(nElm);
    if(connecEdges.size() == 0)
      _connecEdges.resize(nElm * nEdg);

    std::sort(connecNodes.begin(), connecNodes.end());
    _nNod = std::unique(connecNodes.begin(), connecNodes.end()) - connecNodes.begin();

    // Color the elements for partitioning
    colorElements(1);
  };
  ~feCncGeo(){ }

  const std::string &getID() const { return _ID; }
  int getTag() const { return _tag; }
  int getDim() const { return _dim; }

  geometryType getGeometry() const { return _geometry; };
  geometricInterpolant getInterpolant() const { return _interpolant; };

  int getNbNodes() { return _nNod; }
  int getNbNodePerElem() const { return _nNodPerElm; }
  int getNbElm() { return _nElm; }
  int getNbEdgePerElem() { return _nEdg; }

  feSpace *getFeSpace() const { return _space; } // Interpolant geometrique
  feStatus setQuadratureRule(feQuadrature *rule);

  void setMeshPtr(feMesh *mesh) { _mesh = mesh; }

  std::vector<int> getNodeConnectivityCopy() { return _connecNodes; }
  std::vector<int> &getNodeConnectivityRef() { return _connecNodes; }
  std::vector<int> &getEdgeConnectivityRef() { return _connecEdges; }
  std::vector<int> &getElemConnectivityRef() { return _connecElem; }

  int getNodeConnectivity(int iNode) { return _connecNodes[iNode]; }
  int getNodeConnectivity(int numElem, int iNode) const
  {
    return _connecNodes[_nNodPerElm * numElem + iNode];
  }
  int getElementConnectivity(int numElem) { return _connecElem[numElem]; }
  int getEdgeConnectivity(int numElem, int iEdge) { return _connecEdges[_nEdg * numElem + iEdge]; }

  void setNodeConnectivity(int numElem, int iNode, int val)
  {
    _connecNodes[_nNodPerElm * numElem + iNode] = val;
  }
  void setElementConnectivity(int numElem, int val) { _connecElem[numElem] = val; }
  void setEdgeConnectivity(int numElem, int iEdge, int val)
  {
    _connecEdges[_nEdg * numElem + iEdge] = val;
  }

  feStatus computeJacobians();
  // double getJacobianAtQuadNode(int numElem, int iQuadNode){ return _J; }
  const std::vector<double> &getJacobians() const { return _J; }

  // void createPatchNode(int nbElm, int nbNodePerElm, std::vector<int> &cncNodes);
  // void createPatchElm(int nbElm, int nbNodePerElm, std::vector<int> &cncNodes);

  void colorElements(int coloringAlgorithm);

  int getNbColor() { return _nbColor; };
  std::vector<int> &getColorElm() { return _elmToColor; };
  std::vector<int> &getNbElmPerColor() { return _nbElmPerColor; };
  std::vector<std::vector<int> > &getListElmPerColor() { return _listElmPerColor; };
  int getNbElmPerColorI(int i) { return _nbElmPerColor[i]; };
  std::vector<int> &getListElmPerColorI(int i) { return _listElmPerColor[i]; };
  int getElmColored(int iColor, int iElmC) { return _listElmPerColor[iColor][iElmC]; };
};

#endif
