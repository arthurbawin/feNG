#ifndef _FEMESH_
#define _FEMESH_

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <set>

#include "feCncGeo.h"
#include "feSpace.h"
#include "feVertex.h"
#include "feEdge.h"
#include "feElement.h"

class feMesh{

protected:
  std::string _ID;
  int _dim;
  int _nNod;
  int _nEdg;
  int _nTotalElm;
  int _nInteriorElm;
  int _nBoundaryElm;

  std::vector<double> _coord;
  std::vector<Vertex> _vertices;

  int _nCncGeo;
  std::vector<feCncGeo*> _cncGeo;
  std::map<std::string,int> _cncGeoMap;

public:
	feMesh(int nNod = 0, int dim = 0, int nCncGeo = 0, std::string ID = "")
    : _ID(ID), _dim(dim), _nNod(nNod), _nEdg(0), _nCncGeo(nCncGeo)
  {

	};
	virtual ~feMesh() {}

  std::string getID(){ return _ID; }
  int getDim(){ return _dim; }
  int getNbNodes(){ return _nNod; }
  int getNbEdges(){ return _nEdg; }
  int getNbElems(){ return _nTotalElm; }
  int getNbInteriorElems(){ return _nInteriorElm; }
  int getNbBoundaryElems(){ return _nBoundaryElm; }

  std::vector<double> &getCoord(){ return _coord; }
  double getCoordDim(int i, int dim){ return _coord[i*_dim+dim]; }
  double getXCoord(int i){ return _coord[i*_dim]; }
  double getYCoord(int i){ return _coord[i*_dim+1]; }
  double getZCoord(int i){ return _coord[i*_dim+2]; }
  std::vector<double> getCoord(std::string cncGeoID, int numElm);
  std::vector<double> getCoord(int cncGeoTag, int numElm);

  Vertex* getVertex(int iVertex){ return &_vertices[iVertex]; }

  int getNbCncGeo() { return _nCncGeo; }
  std::vector<feCncGeo*> &getCncGeo(){ return _cncGeo; }
  int getCncGeoTag(std::string cncGeoID);
  feCncGeo* getCncGeoByName(std::string cncGeoID);
  feCncGeo* getCncGeoByTag(int cncGeoTag);
  int getNbElm(std::string cncGeoID);
  int getNbElm(int cncGeoTag);
  int getNbNodePerElem(std::string cncGeoID);
  int getNbNodePerElem(int cncGeoTag);
  int getVertex(std::string cncGeoID, int numElem, int numVertex);
  int getVertex(int cncGeoTag, int numElem, int numVertex);
  int getElement(std::string cncGeoID, int numElem);
  int getElement(int cncGeoTag, int numElem);

  // int getEdge(std::string cncGeoID, int numElem, int numEdge);
  // int getEdge(int cncGeoTag, int numElem, int numEdge);

  feSpace* getGeometricSpace(std::string cncGeoID);
  feSpace* getGeometricSpace(int cncGeoTag);

  void printInfo();
};

class feMesh1DP1 : public feMesh{

protected:
  int _nElm;
  double _xA, _xB;
  std::string _bndA_ID, _bndB_ID, _domID;
  int _nElmDomain;
  int _nElmBoundary;
  int _nNodDomain;
  int _nNodBoundary;

public:
  feMesh1DP1(double xA, double xB, int nElm, std::string bndA_ID, std::string bndB_ID, std::string domID);
  virtual ~feMesh1DP1();
};

class feMesh2DP1 : public feMesh{
private:
  // A physical entity (a domain)
  typedef struct physicalEntityStruct {
    std::string name;
    uint32_t dim;
    uint32_t tag;
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
    uint32_t dim;
    uint32_t tag;
    uint32_t tagGmsh;
    uint32_t numPhysicalTags;
    std::vector<int> physicalTags;
    uint32_t gmshType;  // Type d'élément dans gmsh
    std::string cncID;
    int nElm;
    int nNodePerElem;
    int nEdgePerElem;
    std::vector<int> connecElem;
    std::vector<int> connecNodes;
    std::vector<int> connecEdges;
    std::vector<int> isBoundedBy;
    std::vector<int> isBoundaryOf;
  } entity;

public:
  typedef std::map<std::pair<int,int>,std::string> mapType; 

protected:
  int _nPhysicalEntities;
  std::map<std::pair<int,int>,entity> _entities;
  std::vector<physicalEntity> _physicalEntities;

  mapType _physicalEntitiesDescription;

  int _nPoints;
  int _nCurves;
  int _nSurfaces;
  int _nVolumes;

  // std::vector<Vertex> _vertices;
  std::set<Edge, EdgeLessThan> _edges;

  int _nElm;

  bool _isBinary;
  double _gmshVersion;

  // std::vector<feElement*> _elements;
  // std::vector<feElement*> _boundaryElements;

public:
  feMesh2DP1(std::string meshName, bool curved, mapType physicalEntitiesDescription = mapType());
  virtual ~feMesh2DP1();

  int readMsh2(std::istream &input, bool curved, mapType physicalEntitiesDescription);
  int readMsh4(std::istream &input, bool curved);
  int readGmsh(std::string meshName, bool curved, mapType physicalEntitiesDescription = mapType());

  bool isMeshFileBinary(){ return _isBinary; }
  double getGmshVersion(){ return _gmshVersion; }

  mapType getPhysicalEntitiesDescription(){ return _physicalEntitiesDescription; }
};

#endif