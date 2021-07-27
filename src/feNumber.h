#ifndef _FENUMBER_
#define _FENUMBER_

#include <algorithm>

#include "feMesh.h"
#include "feSpace.h"
#include "feCncGeo.h"

#define ESS -1
#define INC -2

class feNumber {
protected:
  std::string _fieldID; // Le nom du champ dont c'est la numerotation

  int _nInc;
  int _nDofs;
  int _nNod;
  int _nElm;
  int _nEdg;

  int _maxDOFperElem;
  int _maxDOFperEdge;

  std::vector<int> _nDOFVertices; // NUMER_ddl_sommet
  std::vector<int> _nDOFElements; // NUMER_ddl_element
  std::vector<int> _nDOFEdges;
  std::vector<int> _numberingVertices; // NUMER_sommet      TODO : reflechir pour l'indexage
  std::vector<int> _numberingElements; // NUMER_element
  std::vector<int> _numberingEdges;

  std::vector<int> _codeDOFVertices; // NUMER_code_ddl_sommet
  std::vector<int> _codeDOFElements; // NUMER_code_ddl_element
  std::vector<int> _codeDOFEdges;

  std::vector<int> _elemToEdge;

public:
  feNumber(feMesh *mesh);
  ~feNumber() {}

  int getNbNodes() { return _nNod; }

  void defDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF);
  void defDDLEdge(feMesh *mesh, std::string cncGeoID, int numElem, int numEdge, int numDOF);
  void defDDLSommet_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem);
  void defDDLEdge_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numEdge);

  int getDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  int getDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF);
  int getDDLEdge(feMesh *mesh, std::string cncGeoID, int numElem, int numEdge, int numDOF);

  void prepareNumbering();
  int numberUnknowns(int globalNum);
  int numberEssential(int globalNum);

  int getVertexNumber(int iVertex) { return _numberingVertices[iVertex]; };

  void printNumberingVertices() {
    printf("nNod = %d\n", _nNod);
    for(size_t i = 0; i < _numberingVertices.size(); ++i) {
      if(_codeDOFVertices[i] == ESS)
        printf("%d \t %s\n", _numberingVertices[i], "ESS");
      else if(_codeDOFVertices[i] == INC)
        printf("%d \t %s\n", _numberingVertices[i], "INC");
      else
        printf("%d \t %d\n", _numberingVertices[i], _codeDOFVertices[i]);
    }
  }

  void printNumberingElements() {
    printf("nElm = %d\n", _nElm);
    printf("maxDOFPerElem = %d\n", _maxDOFperElem);
    for(int i = 0; i < _nElm; ++i) {
      for(int j = 0; j < _maxDOFperElem; ++j) {
        if(_codeDOFElements[i] == ESS)
          printf("elm %d - dof %d : %d \t %s\n", i, j, _numberingElements[_maxDOFperElem * i + j],
                 "ESS");
        else if(_codeDOFElements[i] == INC)
          printf("elm %d - dof %d : %d \t %s\n", i, j, _numberingElements[_maxDOFperElem * i + j],
                 "INC");
        else
          printf("elm %d - dof %d : %d \t %d\n", i, j, _numberingElements[_maxDOFperElem * i + j],
                 _codeDOFElements[i]);
      }
    }
  }

  void printNumberingEdges() {
    printf("nEdg = %d\n", _nEdg);
    printf("maxDOFperEdge = %d\n", _maxDOFperEdge);
    for(int i = 0; i < _nEdg; ++i) {
      for(int j = 0; j < _maxDOFperEdge; ++j) {
        if(_codeDOFEdges[i] == ESS)
          printf("%d \t %s\n", _numberingEdges[_maxDOFperEdge * i + j], "ESS");
        else if(_codeDOFEdges[i] == INC)
          printf("%d \t %s\n", _numberingEdges[_maxDOFperEdge * i + j], "INC");
        else
          printf("%d \t %d\n", _numberingEdges[_maxDOFperEdge * i + j], _codeDOFEdges[i]);
      }
    }
  }

  void printCodeVertices() {
    for(auto val : _codeDOFVertices) std::cout << val << std::endl;
  }

  void printCodeElements() {
    for(auto val : _codeDOFElements) std::cout << val << std::endl;
  }

  void printCodeEdges() {
    for(auto val : _codeDOFEdges) std::cout << val << std::endl;
  }
};

class feMetaNumber {
protected:
  int _nInc;
  int _nDofs;
  int _nFields;
  std::vector<std::string> _fieldIDs;
  std::map<std::string, feNumber *> _numberings;

public:
  feMetaNumber(feMesh *mesh, const std::vector<feSpace *> &space,
               const std::vector<feSpace *> &essBC);
  ~feMetaNumber();

  int getNbUnknowns() { return _nInc; }
  int getNbDOFs() { return _nDofs; }
  int getNbFields() { return _nFields; }
  std::string getFieldID(int iField) { return _fieldIDs[iField]; }

  feNumber *getNumbering(std::string fieldID) {
    return _numberings[fieldID];
  } // TODO : check bounds (-:
  feNumber *getNumbering(int fieldTag) {
    return _numberings[_fieldIDs[fieldTag]];
  } // TODO : check bounds (-:

  void printFields() {
    for(auto s : _fieldIDs) std::cout << s << std::endl;
  }

  void printNumberings() {
    for(int i = 0; i < _nFields; ++i) {
      // for(auto const& secondIsNumber : _numberings){
      std::cout << "Field " << _fieldIDs[i] << " - vertices :" << std::endl;
      _numberings[_fieldIDs[i]]->printNumberingVertices();
      std::cout << "Field " << _fieldIDs[i] << " - elements :" << std::endl;
      _numberings[_fieldIDs[i]]->printNumberingElements();
      std::cout << "Field " << _fieldIDs[i] << " - edges :" << std::endl;
      _numberings[_fieldIDs[i]]->printNumberingEdges();
    }
  }

  void printCodes() { // Attention : ordre dépend du nom des variables (ordre du map.first)
    for(auto const &secondIsNumber : _numberings) {
      std::cout << "Field " << secondIsNumber.first << " - vertices :" << std::endl;
      secondIsNumber.second->printCodeVertices();
      std::cout << "Field " << secondIsNumber.first << " - elements :" << std::endl;
      secondIsNumber.second->printCodeElements();
      std::cout << "Field " << secondIsNumber.first << " - edges :" << std::endl;
      secondIsNumber.second->printCodeEdges();
    }
  }
};

#endif