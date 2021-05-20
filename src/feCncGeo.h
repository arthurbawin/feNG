#ifndef _FECNCGEO_
#define _FECNCGEO_

class feSpace;

class feCncGeo{

protected:
  std::string _ID;
  std::string _forme; // ?
  int _nNod;                   // Nombre de noeuds par element
  int _nElm;                   // Nombre d'elements dans cette connectivite
  int _nEdg;
  std::vector<int> _connecNodes; // Connecvitite locale
  std::vector<int> _connecElem;  //_nmCnc;  // Connecvitite globale 
  std::vector<int> _connecEdges;
  std::vector<int> _connecFaces;

  feSpace *_space;

public:
  feCncGeo(int nNod, int nElm, int nEdg, std::string ID, std::string forme, feSpace *space, std::vector<int> connecNodes, 
    std::vector<int> connecElem  = std::vector<int>(), std::vector<int> connecEdges = std::vector<int>(),
    std::vector<int> connecFaces = std::vector<int>()) :  _ID(ID), _forme(forme), _nNod(nNod), _nElm(nElm), _nEdg(nEdg),
    _connecNodes(connecNodes),  _connecElem(connecElem), _connecEdges(connecEdges),  _connecFaces(connecFaces),  _space(space)
  {
    if(connecElem.size() == 0) _connecElem.resize(nElm);
  };
	~feCncGeo(){
    // if(_space != nullptr)
    //   delete _space;
  }

  std::string getID(){ return _ID; }
  std::string getForme(){ return _forme; }
  int getNbNodePerElem(){ return _nNod; }
  int getNbElm()        { return _nElm; }
  int getNbEdgePerElem(){ return _nEdg; }
  feSpace* getFeSpace(){ return _space; } // Interpolant geometrique

  int getLocalConnectivity(int numElem, int iNode){ return _connecNodes[_nNod*numElem+iNode]; }
  int getNodeConnectivity(int numElem, int iNode){ return _connecNodes[_nNod*numElem+iNode]; }
  int getGlobalConnectivity(int numElem)          { return _connecElem[numElem]; }
  int getEdgeConnectivity(int numElem, int iEdge) { return _connecEdges[_nEdg*numElem+iEdge]; }

  void setLocalConnectivity(int numElem, int iNode, int val){ _connecNodes[_nNod*numElem+iNode] = val; }
  void setNodeConnectivity(int numElem, int iNode, int val){ _connecNodes[_nNod*numElem+iNode] = val; }
  void setGlobalConnectivity(int numElem, int val)          { _connecElem[numElem] = val; }
  void setEdgeConnectivity(int numElem, int iEdge, int val) { _connecEdges[_nEdg*numElem+iEdge] = val; }
};

#endif