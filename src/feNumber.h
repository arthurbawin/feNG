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
public:
  feNumber(feMesh *mesh);
  ~feNumber() {}

  int getNbNodes(){ return _nNod; }

  int getDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  int getDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLSommet(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement(feMesh *mesh, std::string cncGeoID, int numElem, int numDOF);
  void defDDLEdge(feMesh *mesh, std::string cncGeoID, int numElem, int numEdge, int numDOF);
  void defDDLSommet_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numVertex);
  void defDDLElement_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem);
  void defDDLEdge_essentialBC(feMesh *mesh, std::string cncGeoID, int numElem, int numEdge);

  void prepareNumbering();
  int numberUnknowns(int globalNum);
  int numberEssential(int globalNum);

  int getVertexNumber(int iVertex){ return _numberingVertices[iVertex]; };

  void printNumberingVertices(){
    for(auto val : _numberingVertices)
      std::cout<<val<<std::endl;
  }

  void printNumberingElements(){
    for(auto val : _numberingElements)
      std::cout<<val<<std::endl;
  }

  void printNumberingEdges(){
    for(auto val : _numberingEdges)
      std::cout<<val<<std::endl;
  }

};

class feMetaNumber{

protected:
	int _nInc;
	int _nDofs;
	int _nFields;
  std::vector<std::string> _fieldIDs;
  std::map<std::string,feNumber*> _numberings;

public:
	feMetaNumber(feMesh *mesh, const std::vector<feSpace*> &space, const std::vector<feSpace*> &essBC);
	~feMetaNumber();

  int getNbUnknowns(){ return _nInc; }
  int getNbDOFs(){ return _nDofs; }
  int getNbFields(){ return _nFields; }
  std::string getFieldID(int iField){ return _fieldIDs[iField]; }

  feNumber* getNumbering(std::string fieldID){ return _numberings[fieldID]; } // TODO : check bounds (-:
  feNumber* getNumbering(int fieldTag){ return _numberings[_fieldIDs[fieldTag]]; } // TODO : check bounds (-:

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
      std::cout<<"Field "<<secondIsNumber.first<<" - edges :"<<std::endl;
      secondIsNumber.second->printNumberingEdges();
    }
  }
};

#endif