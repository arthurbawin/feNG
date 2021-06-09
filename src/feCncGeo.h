#ifndef _FECNCGEO_
#define _FECNCGEO_

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class feSpace;
class feMesh;

class feCncGeo{

protected:
  std::string _ID;
  int _tag;
  int _dim;
  std::string _forme;
  int _nNod;
  int _nNodPerElm;             // Nombre de noeuds par element
  int _nElm;                   // Nombre d'elements dans cette connectivite
  int _nEdg;
  std::vector<int> _connecNodes; // Connectivite des noeuds
  std::vector<int> _connecElem;  // Connectivite des elements
  std::vector<int> _connecEdges; // Connectivite des aretes
  std::vector<int> _connecFaces; // Connectivite des faces

  feSpace *_space;
  feMesh *_mesh; // Ou alors donner le vecteur de Vertex de la cnc

  std::vector<double> _J; // Jacobiens

public:
  feCncGeo(int tag, int dim, int nNod, int nElm, int nEdg, std::string ID, std::string forme, feSpace *space, 
    std::vector<int> connecNodes, 
    std::vector<int> connecElem  = std::vector<int>(), 
    std::vector<int> connecEdges = std::vector<int>(), 
    std::vector<int> connecFaces = std::vector<int>())
    : _ID(ID), _tag(tag), _dim(dim), _forme(forme), _nNodPerElm(nNod), _nElm(nElm), _nEdg(nEdg),
    _connecNodes(connecNodes),  _connecElem(connecElem), _connecEdges(connecEdges),  _connecFaces(connecFaces),  _space(space)
  {
    if(connecElem.size() == 0) _connecElem.resize(nElm);
    std::sort(connecNodes.begin(), connecNodes.end());
    _nNod = std::unique(connecNodes.begin(), connecNodes.end()) - connecNodes.begin();
  };
	~feCncGeo(){
    // if(_space != nullptr)
    //   delete _space;
  }

  std::string getID()   { return _ID;         }
  int getDim()          { return _dim;        }
  std::string getForme(){ return _forme;      }
  int getNbNodes()      { return _nNod;       }
  int getNbNodePerElem(){ return _nNodPerElm; }
  int getNbElm()        { return _nElm;       }
  int getNbEdgePerElem(){ return _nEdg;       }
  feSpace* getFeSpace() { return _space;      } // Interpolant geometrique

  void setMeshPtr(feMesh *mesh){ _mesh = mesh; }

  std::vector<int> getNodeConnectivityCopy(){ return _connecNodes; }
  std::vector<int> &getNodeConnectivityRef(){ return _connecNodes; }
  std::vector<int> &getElemConnectivityRef(){ return _connecElem;  }

  int getLocalConnectivity(int numElem, int iNode){ return _connecNodes[_nNodPerElm*numElem+iNode]; }
  int getNodeConnectivity(int numElem, int iNode) { return _connecNodes[_nNodPerElm*numElem+iNode]; }
  int getGlobalConnectivity(int numElem)          { return _connecElem[numElem]; }
  int getElementConnectivity(int numElem)         { return _connecElem[numElem]; }
  int getEdgeConnectivity(int numElem, int iEdge) { return _connecEdges[_nEdg*numElem+iEdge]; }

  void setLocalConnectivity(int numElem, int iNode, int val){ _connecNodes[_nNodPerElm*numElem+iNode] = val; }
  void setNodeConnectivity(int numElem, int iNode, int val) { _connecNodes[_nNodPerElm*numElem+iNode] = val; }
  void setGlobalConnectivity(int numElem, int val)          { _connecElem[numElem] = val; }
  void setElementConnectivity(int numElem, int val)         { _connecElem[numElem] = val; }
  void setEdgeConnectivity(int numElem, int iEdge, int val) { _connecEdges[_nEdg*numElem+iEdge] = val; }

  void computeJacobians();
  // double getJacobianAtQuadNode(int numElem, int iQuadNode){ return _J; }
  std::vector<double> &getJacobians(){ return _J; }
};

#endif