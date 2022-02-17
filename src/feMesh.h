#ifndef _FEMESH_
#define _FEMESH_

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <set>

#include "feMessage.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feVertex.h"
#include "feEdge.h"
#include "feElement.h"
#include "feTriangle.h"
#include "rtree.h"

class feMesh
{
protected:
  std::string _ID;
  int _dim;
  int _nNod;
  int _nNodesWithNoPhysical;
  int _nEdg;
  int _nTotalElm;
  int _nInteriorElm;
  int _nBoundaryElm;

  std::vector<double> _coord;
  std::vector<Vertex> _vertices;
  std::map<int, int> _verticesMap;

  int _nCncGeo;
  std::vector<feCncGeo *> _cncGeo;

  std::map<std::string, int> _cncGeoMap;

public:
  feMesh(int nNod = 0, int dim = 0, int nCncGeo = 0, std::string ID = "")
    : _ID(ID), _dim(dim), _nNod(nNod), _nEdg(0), _nCncGeo(nCncGeo){};
  virtual ~feMesh() {}

  std::string getID() { return _ID; }
  int getDim() { return _dim; }
  int getNbNodes() { return _nNod; }
  int getNbNodesWithoutPhysical() { return _nNodesWithNoPhysical; }
  int getNbEdges() { return _nEdg; }
  int getNbElems() { return _nTotalElm; }
  int getNbInteriorElems() { return _nInteriorElm; }
  int getNbBoundaryElems() { return _nBoundaryElm; }

  std::vector<double> &getCoord() { return _coord; }
  double getCoordDim(int i, int dim) { return _coord[i * _dim + dim]; }
  double getXCoord(int i) { return _coord[i * _dim]; }
  double getYCoord(int i) { return _coord[i * _dim + 1]; }
  double getZCoord(int i) { return _coord[i * _dim + 2]; }
  std::vector<double> getCoord(std::string cncGeoID, int numElm);
  std::vector<double> getCoord(int cncGeoTag, int numElm);

  int getVertexSequentialTagFromGmshTag(int gmshNodeTag) { return _verticesMap[gmshNodeTag]; }
  Vertex *getVertex(int iVertex) { return &_vertices[iVertex]; }
  Vertex *getVertexFromGmshNodeTag(int gmshNodeTag)
  {
    return &_vertices[_verticesMap[gmshNodeTag]];
  }
  std::vector<Vertex> &getVertices() { return _vertices; }

  std::set<Edge, EdgeLessThan> _edges;
  std::vector<Triangle *> _elements; // Hardcoded triangles for now

  int getNbCncGeo() { return _nCncGeo; }
  std::vector<feCncGeo *> &getCncGeo() { return _cncGeo; }
  int getCncGeoTag(std::string cncGeoID);
  feCncGeo *getCncGeoByName(std::string cncGeoID);
  feCncGeo *getCncGeoByTag(int cncGeoTag);
  int getNbElm(std::string cncGeoID);
  int getNbElm(int cncGeoTag);
  int getNbNodePerElem(std::string cncGeoID);
  int getNbNodePerElem(int cncGeoTag);
  int getVertex(std::string cncGeoID, int numElem, int numVertex);
  int getVertex(int cncGeoTag, int numElem, int numVertex);
  int getElement(std::string cncGeoID, int numElem);
  int getElement(int cncGeoTag, int numElem);
  int getEdge(std::string cncGeoID, int numElem, int numEdge);
  int getEdge(int cncGeoTag, int numElem, int numEdge);
  feSpace *getGeometricSpace(std::string cncGeoID);
  feSpace *getGeometricSpace(int cncGeoTag);

  virtual bool locateVertex(std::vector<double> &x, int &iElm, std::vector<double> &u,
                            double tol = 1e-5) = 0;

  void printInfo(bool printConnectivities = true);

  std::vector<int> color(int cncGeoTag = 0); //i=0 => borders; i=1 =>element
};

class feMesh1DP1 : public feMesh
{
protected:
  int _nElm;
  double _xA, _xB;
  std::string _bndA_ID, _bndB_ID, _domID;
  int _nElmDomain;
  int _nElmBoundary;
  int _nNodDomain;
  int _nNodBoundary;

public:
  feMesh1DP1(double xA, double xB, int nElm, std::string bndA_ID, std::string bndB_ID,
             std::string domID);
  virtual ~feMesh1DP1();

  virtual bool locateVertex(std::vector<double> &x, int &iElm, std::vector<double> &u,
                            double tol = 1e-5){};
};

class feMesh0DP0 : public feMesh
{
protected:
  int _nElm;
  double _xA;
  std::string _domID;
  int _nElmDomain;
  int _nElmBoundary;
  int _nNodDomain;
  int _nNodBoundary;

public:
  feMesh0DP0(double xA, int nElm, std::string domID);
  virtual ~feMesh0DP0();

  virtual bool locateVertex(std::vector<double> &x, int &iElm, std::vector<double> &u,
                            double tol = 1e-5){};
};

class feMetaNumber;
class feSolutionContainer;

class feMesh2DP1 : public feMesh
{
private:
  // A physical entity (a domain)
  typedef struct physicalEntityStruct {
    std::string name;
    int dim;
    int tag;
    std::vector<int> listEntities; // Tag des entities
    std::string cncID;
    feSpace *geoSpace;
    int nElm;
    int nNodePerElem;
    int nEdgePerElem;
    std::vector<int> connecElem;
    std::vector<int> connecNodes;
    std::vector<int> connecEdges;
    std::vector<int> isBoundedBy;
    std::vector<int> isBoundaryOf;
  } physicalEntity;
  // A geometric entity (building block to make physical entities)
  typedef struct entityStruct {
    int dim;
    int tag;
    int tagGmsh;
    int numPhysicalTags;
    std::vector<int> physicalTags;
    int gmshType; // Type d'élément dans gmsh
    std::string cncID;
    int nElm;
    int nNodePerElem = 0;
    int nEdgePerElem;
    std::vector<int> connecElem;
    std::vector<int> connecNodes;
    std::vector<int> connecEdges;
    std::vector<int> isBoundedBy;
    std::vector<int> isBoundaryOf;
  } entity;

  RTree<int, double, 3> _rtree;

public:
  // FIXME : Choose more explicit name, for example "physicalEntitiesDescription"
  typedef std::map<std::pair<int, int>, std::string> mapType;

protected:
  int _nPhysicalEntities;
  std::map<std::pair<int, int>, entity> _entities;
  std::map<std::pair<int, int>, physicalEntity> _physicalEntities;

  mapType _physicalEntitiesDescription;

  std::vector<int> _globalCurvesNodeConnectivity;
  std::vector<int> _globalSurfacesNodeConnectivity;
  std::vector<int> _sequentialNodeToGmshNode;

  int _nPoints;
  int _nCurves;
  int _nSurfaces;
  int _nVolumes;

  // std::vector<Vertex> _vertices;
  // std::set<Edge, EdgeLessThan> _edges;

  int _nElm;

  bool _isBinary;
  double _gmshVersion;

  // std::vector<Triangle *> _elements; // Hardcoded triangles for now
  // std::vector<feElement*> _elements;
  // std::vector<feElement*> _boundaryElements;

public:
  // Empty constructor : the mesh must be created afterwards with feCheck(readGmsh(...))
  feMesh2DP1(){};
  // Constructor calling readGmsh : will probably be removed
  feMesh2DP1(std::string meshName, bool curved = false, bool reversed = false,
             mapType physicalEntitiesDescription = mapType());
  ~feMesh2DP1();

  feStatus readMsh2(std::istream &input, bool curved, bool reversed, mapType physicalEntitiesDescription);
  feStatus readMsh4(std::istream &input, bool curved, bool reversed, mapType physicalEntitiesDescription);
  feStatus readGmsh(std::string meshName, bool curved = false, bool reversed = false,
                    mapType physicalEntitiesDescription = mapType());

  bool isMeshFileBinary() { return _isBinary; }
  double getGmshVersion() { return _gmshVersion; }

  mapType getPhysicalEntitiesDescription() { return _physicalEntitiesDescription; }

  int getGmshNodeTag(int iVertex) { return _sequentialNodeToGmshNode[iVertex]; }

  bool locateVertex(std::vector<double> &x, int &iElm, std::vector<double> &u, double tol = 1e-5);

  void transfer(feMesh2DP1 *otherMesh, feMetaNumber *myMN, feMetaNumber *otherMN,
                feSolutionContainer *solutionContainer, const std::vector<feSpace *> &mySpaces,
                const std::vector<feSpace *> &mySpacesEssBC,
                const std::vector<feSpace *> &otherSpaces);


 
};

#endif
