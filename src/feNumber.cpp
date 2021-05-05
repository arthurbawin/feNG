#include "feNumber.h"

feNumber::feNumber(feMesh *mesh) : _nNod(mesh->getNbNodes()){
  _nElm = 0;
  for(int i = 0; i < mesh->getCncGeo().size(); ++i)
    _nElm += mesh->getCncGeo()[i]->getNbElm();

  _nDOFVertices.resize(_nNod);
  _nDOFElements.resize(_nElm);
  _codeDOFVertices.resize(_nNod);
  _codeDOFElements.resize(_nElm);
};

int feNumber::getDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex){
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  return _numberingVertices[vert];
}

int feNumber::getDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF){
  int elem = mesh->getElement(cncGeoID, numElem);
  return _numberingElements[_maxDOFperElem * numElem + numDOF];
}

void feNumber::defDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex){
   int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
   _nDOFVertices[vert] = 1;
   _codeDOFVertices[vert] = INC;
}

void feNumber::defDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF){
   int elem = mesh->getElement(cncGeoID, numElem);
   _nDOFElements[elem] = numDOF;
   _codeDOFElements[elem] = INC;
}

void feNumber::defDDLSommet_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex){
   int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
   _codeDOFVertices[vert] = ESS;
}

void feNumber::defDDLElement_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem){
   int elem = mesh->getElement(cncGeoID, numElem);
   _codeDOFElements[elem] = ESS;
}

// TODO : defDDLElement pour les P2 et +

void feNumber::prepareNumbering(){
  _maxDOFperElem = *std::max_element(_nDOFElements.begin(), _nDOFElements.end());
  _numberingVertices.resize(_nNod);
  _numberingElements.resize(_nElm * _maxDOFperElem);

  for(int iVertex = 0; iVertex < _nNod; ++iVertex){
    if(_nDOFVertices[iVertex] == 1)
      _numberingVertices[iVertex] = _codeDOFVertices[iVertex];
  }
  for(int iElm = 0; iElm < _nElm; ++iElm){
    for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
      _numberingElements[_maxDOFperElem * iElm + iDOF] = _codeDOFElements[iElm];
  }
}

int feNumber::numberUnknowns(int globalNum){
  _nInc = 0;
  for(int i = 0; i < _nNod; ++i){
    if(_numberingVertices[i] == INC){
      ++_nInc;
      _numberingVertices[i] = globalNum++;
    }
  }
  for(int i = 0; i < _nElm; ++i){
    for(int j = 0; j < _nDOFElements[i]; ++j){
      if(_numberingElements[_maxDOFperElem * i + j] == INC){
        ++_nInc;
        _numberingElements[_maxDOFperElem * i + j] = globalNum++;
      }
    }
  }
  return globalNum;
}

int feNumber::numberEssential(int globalNum){
  _nDofs = _nInc;
  for(int i = 0; i < _nNod; ++i){
    if(_numberingVertices[i] == ESS){
      ++_nDofs;
      _numberingVertices[i] = globalNum++;
    }
  }
  for(int i = 0; i < _nElm; ++i){
    for(int j = 0; j < _nDOFElements[i]; ++j){
      if(_numberingElements[_maxDOFperElem * i + j] == ESS){
        ++_nDofs;
        _numberingElements[_maxDOFperElem * i + j] = globalNum++;
      }
    }
  }
  return globalNum;
}

feMetaNumber::feMetaNumber(feMesh *mesh, const std::vector<feSpace*> &space, const std::vector<feSpace*> &essBC){
  // AJOUTECHAMP
  for(feSpace *fS : space){
    if(std::find(_fieldIDs.begin(), _fieldIDs.end(), fS->getFieldID()) == _fieldIDs.end())
      _fieldIDs.push_back(fS->getFieldID());
  }
  _nFields = _fieldIDs.size();
  // VRFCHAMP
  for(feSpace *fSBC : essBC){
    bool err = true;
    for(feSpace *fS : space)
      if(fSBC->getFieldID().compare(fS->getFieldID()) == 0 && fSBC->getCncGeoID().compare(fS->getCncGeoID()) == 0)
        err = false;
    if(err) std::cout<<"ATTENTION : Condition essentielle imposée sur un champ ou une connectivité non présent(e)."<<std::endl;
  }
  // Une numerotation pour chaque champ
  for(int i = 0; i < _nFields; ++i)
    _numberings[_fieldIDs[i]] = new feNumber(mesh);

  // for (auto const& x : _numberings)
  //   std::cout << x.first << ':'  << x.second  << std::endl;

  // INITIALISE_FENUMER
  for(feSpace *fS : space){
    fS->initializeNumberingUnknowns(_numberings[fS->getFieldID()]);
  }
  // INITIALISE_FENUMER_CLMESS
  for(feSpace *fSBC : essBC){
    fSBC->initializeNumberingEssential(_numberings[fSBC->getFieldID()]);
  }
  // PREPARER_NUMEROTATION
  // for(std::map< >::iterator s : _numberings ){
  //   string fieldID = s.first;
  //   feNumber* number = s.second;
  // }
  for(auto const& secondIsNumber : _numberings)
    secondIsNumber.second->prepareNumbering();
  // NUMEROTER
  int globalNum = 0;
  for(auto const& secondIsNumber : _numberings)
    globalNum = secondIsNumber.second->numberUnknowns(globalNum);
  _nInc = globalNum;
  for(auto const& secondIsNumber : _numberings)
    globalNum = secondIsNumber.second->numberEssential(globalNum);
  _nDofs = globalNum;
};

feMetaNumber::~feMetaNumber() {
  // Delete all numberings
  for(std::map<std::string, feNumber*>::iterator it = _numberings.begin(); it != _numberings.end(); ++it){
    delete it->second; // delete le feNumber associé à it
  }
}