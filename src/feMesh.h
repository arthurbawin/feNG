#ifndef _FEMESH_
#define _FEMESH_

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "feCncGeo.h"
#include "feSpace.h"

class feMesh{

protected:
  std::string _ID;
  int _dim;
  int _nNod;
  int _nTotalElm;

  std::vector<double> _coord;

  int _nCncGeo;
  std::vector<feCncGeo*> _cncGeo;
  std::map<std::string,int> _cncGeoMap;

public:
	feMesh(int nNod, int dim, int nCncGeo, std::string ID)
    : _ID(ID), _dim(dim), _nNod(nNod), _nCncGeo(nCncGeo)
  {

	};
	virtual ~feMesh() {}

  std::string getID(){ return _ID; }
  int getDim(){ return _dim; }
  int getNbNodes() { return _nNod; }
  int getNbElems(){ return _nTotalElm; }

  std::vector<double> &getCoord(){ return _coord; }
  double getCoordDim(int i, int dim){ return _coord[i*_dim+dim]; }
  double getXCoord(int i){ return _coord[i*_dim]; }
  double getYCoord(int i){ return _coord[i*_dim+1]; }
  double getZCoord(int i){ return _coord[i*_dim+2]; }
  std::vector<double> getCoord(std::string cncGeoID, int numElm);
  std::vector<double> getCoord(int cncGeoTag, int numElm);

  int getNbCncGeo() { return _nCncGeo; }
  std::vector<feCncGeo*> getCncGeo(){ return _cncGeo; }
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

  feSpace* getGeometricSpace(std::string cncGeoID);
  feSpace* getGeometricSpace(int cncGeoTag);
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
  virtual ~feMesh1DP1() {
    for(feCncGeo *cnc : _cncGeo){
      delete cnc->getFeSpace();
      delete cnc;
    }
  }

  void printInfo(){
    std::cout<<"Nombre total d'elements : "<<_nTotalElm<<std::endl;
    std::cout<<"Nombre de connectivites : "<<_nCncGeo<<std::endl;
    for(feCncGeo *cnc : _cncGeo){
      std::cout<< "CncGeo "<<cnc->getID()<<std::endl;
      for(int i = 0; i < cnc->getNbElm(); ++i)
        std::cout<< cnc->getGlobalConnectivity(i) << std::endl;
    }
    for(feCncGeo *cnc : _cncGeo){
      std::cout<< "CncGeo "<<cnc->getID()<<std::endl;
      for(int i = 0; i < cnc->getNbElm(); ++i)
        std::cout<< cnc->getGlobalConnectivity(i) << std::endl;
    }
  }

};

#endif