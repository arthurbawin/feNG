#ifndef _FECNCGEO_
#define _FECNCGEO_

class feSpace;

class feCncGeo{

protected:
  std::string _ID;
  std::string _forme; // ?
  int _nNod;                   // Nombre de noeuds par element
  int _nElm;                   // Nombre d'elements dans cette connectivite
  std::vector<int> _connecLoc; // Connecvitite locale
  std::vector<int> _connecGlo; //_nmCnc;  // Connecvitite globale 

  feSpace *_space;

public:
	feCncGeo(int nNod, int nElm, std::vector<int> connec, std::string ID, std::string forme, feSpace *space)
    : _nNod(nNod), _nElm(nElm), _connecLoc(connec), _ID(ID), _forme(forme), _space(space)
  {
    _connecGlo.resize(nElm);
	};
	~feCncGeo() {}

  std::string getID(){ return _ID; }
  int getNbNodePerElem(){ return _nNod; }
  int getNbElm(){ return _nElm; }
  feSpace* getFeSpace(){ return _space; } // Interpolant geometrique

  int getLocalConnectivity(int numElem, int iNode){ return _connecLoc[_nNod*numElem+iNode]; }
  void setLocalConnectivity(int numElem, int iNode, int val){ _connecLoc[_nNod*numElem+iNode] = val; }
  int getGlobalConnectivity(int numElem){ return _connecGlo[numElem]; }
  void setGlobalConnectivity(int numElem, int val){ _connecGlo[numElem] = val; }
};

#endif