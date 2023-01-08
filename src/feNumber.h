#ifndef _FENUMBER_
#define _FENUMBER_

#include <algorithm>

#include "feMesh.h"
#include "feSpace.h"
#include "feCncGeo.h"

typedef enum {DOF_ESSENTIAL = -1, DOF_UNKNOWN = -2, DOF_NOT_ASSIGNED = -3} dofType;

class feNumber
{
protected:
  std::string _fieldID; // Le nom du champ dont c'est la numerotation

  int _nInc;
  int _nDofs;
  int _nNod;
  int _nElm;
  int _nEdg;

  int _maxDOFperVertex;
  int _maxDOFperElem;
  int _maxDOFperEdge;

  std::vector<int> _nDOFVertices; // NUMER_ddl_sommet
  std::vector<int> _nDOFElements; // NUMER_ddl_element
  std::vector<int> _nDOFEdges;
  std::vector<int> _numberingVertices; // NUMER_sommet
  std::vector<int> _numberingElements; // NUMER_element
  std::vector<int> _numberingEdges;

  std::vector<int> _codeDOFVertices; // NUMER_code_ddl_sommet
  std::vector<int> _codeDOFElements; // NUMER_code_ddl_element
  std::vector<int> _codeDOFEdges;

  std::vector<int> _elemToEdge;

  // The index of all the degrees of freedom associated to this field.
  // Elements, faces (edges) and vertices are mixed.
  std::vector<int> _allEssentialDOF;
  std::vector<int> _allUnknownDOF;

public:
  feNumber(feMesh *mesh);
  ~feNumber() {};

  int getNbNodes() { return _nNod; }
  int getNbDOFs() { return _nDofs; }

  // Assign 1 DOF per vertex by default for continuous elements (i.e. not for DG)
  void defDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF = 1);
  void defDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF);
  void defDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge, int numDOF);
  void defDDLSommet_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                int numVertex);
  void defDDLElement_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem);
  void defDDLEdge_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge);

  int getDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF = 0);
  int getDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF);
  int getDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge, int numDOF);

  void prepareNumbering();
  int numberUnknowns(int globalNum);
  int numberEssential(int globalNum);

  std::vector<int> &getEssentialDOF() { return _allEssentialDOF; };
  std::vector<int> &getUnknownDOF() { return _allUnknownDOF; };

  int getDOFNumberAtVertex(int iVertex) { return _numberingVertices[iVertex]; };
  int getDOFNumberAtVertex(int iVertex, int numDOF) { return _numberingVertices[_maxDOFperVertex * iVertex + numDOF]; };
  int getDOFCodeAtVertex(int iVertex) { return _codeDOFVertices[iVertex]; };

  void exportNumberingVertices(feMesh *mesh, FILE *file);

  void compactFieldDOF();

public:
  void printNumberingVertices();
  void printNumberingElements();
  void printNumberingEdges();
  void printCodeVertices();
  void printCodeElements();
  void printCodeEdges();
};

class feMetaNumber
{
protected:
  int _nInc;
  int _nDofs;
  int _nFields;
  std::vector<std::string> _fieldIDs;
  std::map<std::string, feNumber*> _numberings;

public:
  feMetaNumber(feMesh *mesh, const std::vector<feSpace *> &space,
               const std::vector<feSpace *> &essBC);
  ~feMetaNumber();

  int getNbUnknowns() { return _nInc; }
  int getNbDOFs() { return _nDofs; }
  int getNbDOFs(std::string fieldID) { return _numberings[fieldID]->getNbDOFs(); }
  int getNbFields() { return _nFields; }
  std::string getFieldID(int iField) { return _fieldIDs[iField]; }

  feNumber *getNumbering(std::string fieldID)
  {
    return _numberings[fieldID];
  } // TODO : check bounds (-:
  feNumber *getNumbering(int fieldTag)
  {
    return _numberings[_fieldIDs[fieldTag]];
  } // TODO : check bounds (-:

  std::vector<int> &getEssentialDOF(std::string fieldID) { return _numberings[fieldID]->getEssentialDOF(); };
  std::vector<int> &getEssentialDOF(int fieldTag) { return _numberings[_fieldIDs[fieldTag]]->getEssentialDOF(); };
  std::vector<int> &getUnknownDOF(std::string fieldID) { return _numberings[fieldID]->getUnknownDOF(); };
  std::vector<int> &getUnknownDOF(int fieldTag) { return _numberings[_fieldIDs[fieldTag]]->getUnknownDOF(); };

  feStatus exportNumberingVertices(feMesh *mesh, std::string fileName);

  void printFields();
  void printNumberings();
  void printCodes();
};

#endif