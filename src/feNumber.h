#ifndef _FENUMBER_
#define _FENUMBER_

#include <algorithm>

#include "feMesh.h"
#include "feSpace.h"
#include "feCncGeo.h"

#define ESS -1
#define INC -2

class feNumber{
protected:
  std::string _fieldID; // Le nom du champ dont c'est la numerotation

  int _nInc;
  int _nDofs;
  int _nNod;
  int _nElm;

  int _maxDOFperElem;

  std::vector<int> _nDOFVertices; // NUMER_ddl_sommet
  std::vector<int> _nDOFElements; // NUMER_ddl_element
  std::vector<int> _numberingVertices; // NUMER_sommet      TODO : reflechir pour l'indexage
  std::vector<int> _numberingElements; // NUMER_element

  std::vector<int> _codeDOFVertices; // NUMER_code_ddl_sommet
  std::vector<int> _codeDOFElements; // NUMER_code_ddl_element
public:
  feNumber(feMesh *mesh);
  ~feNumber() {}

  int getDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  int getDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF);
  void defDDLSommet_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem);

  void prepareNumbering();
  int numberUnknowns(int globalNum);
  int numberEssential(int globalNum);

  void printNumberingVertices(){
    for(auto val : _numberingVertices)
      std::cout<<val<<std::endl;
  }

  void printNumberingElements(){
    for(auto val : _numberingElements)
      std::cout<<val<<std::endl;
  }

};

class feMetaNumber{

protected:
	int _nInc;
	int _nDofs;
	int _nFields;
  std::vector<std::string> _fieldIDs;
  std::map<std::string,feNumber*> _numberings; // Destructeur par defaut

public:
	feMetaNumber(feMesh *mesh, const std::vector<feSpace*> &space, const std::vector<feSpace*> &essBC);
	~feMetaNumber();

  int getNbUnknowns(){ return _nInc; }
  int getNbDOFs(){ return _nDofs; }

  feNumber* getNumbering(std::string fieldID){ return _numberings[fieldID]; } // TODO : check bounds (-:

  void printFields(){
    for(auto s : _fieldIDs)
      std::cout<<s<<std::endl;
  }

  void printNumberings(){
    for(auto const& secondIsNumber : _numberings){
      std::cout<<"Field "<<secondIsNumber.first<<" - vertices :"<<std::endl;
      secondIsNumber.second->printNumberingVertices();
      std::cout<<"Field "<<secondIsNumber.first<<" - elements :"<<std::endl;
      secondIsNumber.second->printNumberingElements();
    }
  }
};

#endif