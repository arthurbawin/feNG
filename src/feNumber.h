#ifndef _FENUMBER_
#define _FENUMBER_

#include <algorithm>

#include "feMesh.h"
#include "feSpace.h"
#include "feCncGeo.h"

typedef enum { DOF_ESSENTIAL = -1, DOF_UNKNOWN = -2, DOF_NOT_ASSIGNED = -3 } dofType;

//
// Handles the numbering of the degrees of freedom (DOF) of a single field.
//
class feNumber
{
protected:
  // Field name
  std::string _fieldID;

  // Number of unknown DOF in this field
  int _nInc;
  // Total number of DOF in this field
  int _nDofs;
  // A copy of the number of vertices, elements, edges and faces in the mesh
  int _nNod;
  int _nElm;
  int _nEdg;
  int _nFac;

  // Maximum number of DOF on a vertex, element or edge
  int _maxDOFperVertex;
  int _maxDOFperElem;
  int _maxDOFperEdge;
  int _maxDOFperFace;

  // Acutal number of DOF at each vertex/element/edge/face
  std::vector<int> _nDOFVertices;
  std::vector<int> _nDOFElements;
  std::vector<int> _nDOFEdges;
  std::vector<int> _nDOFFaces;

  // Arrays of size nNod/nElm/nEdg x (max # DOF)
  // Contains the global number of the DOF or DOF_NOT_ASSIGNED is
  // the DOF is not used ("ghost" DOF)
  std::vector<int> _numberingVertices;
  std::vector<int> _numberingElements;
  std::vector<int> _numberingEdges;
  std::vector<int> _numberingFaces;

  // Code (see enum above) of each DOF
  std::vector<int> _codeDOFVertices;
  std::vector<int> _codeDOFElements;
  std::vector<int> _codeDOFEdges;
  std::vector<int> _codeDOFFaces;

  std::vector<int> _elemToEdge;

  // The global number of all essential DOF (elements, faces (edges) and vertices together).
  std::set<int> _allEssentialDOF;
  // The global number of all unknown DOF (elements, faces (edges) and vertices together).
  std::set<int> _allUnknownDOF;

public:
  friend class feMetaNumber;

  // Initialize and allocate vectors. Should not be used directly, but called
  // through the feMetaNumber class.
  feNumber(feMesh *mesh, const std::string &fieldName);
  ~feNumber(){};

  // int getNbNodes() { return _nNod; }
  int getNbDOFs() { return _nDofs; }

  // Return a reference to the set of all UNKNOWN DOF of this field (vert/elem/edge mixed)
  std::set<int> &getUnknownDOF() { return _allUnknownDOF; };
  // Return a reference to the set of all ESSENTIAL DOF of this field (vert/elem/edge mixed)
  std::set<int> &getEssentialDOF() { return _allEssentialDOF; };

  // Set the numDOF-th DOF associated to the numVertex-th vertex on the numElem-th element
  // of geometric connectivity cncGeoID as UNKNOWN. The default number of DOF per vertex
  // for continuous elements (not DG) is 1.
  void setUnknownVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numVertex,
                           int numDOF = 1);
  // Set the numDOF-th DOF associated to the numElem-th element
  // of geometric connectivity cncGeoID as UNKNOWN.
  void setUnknownElementDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numDOF);
  // Set the numDOF-th DOF associated to the numEdge-th edge on the numElem-th element
  // of geometric connectivity cncGeoID as UNKNOWN.
  void setUnknownEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numEdge,
                         int numDOF);
  // Set the numDOF-th DOF associated to the numFace-th face on the numElem-th element
  // of geometric connectivity cncGeoID as UNKNOWN.
  void setUnknownFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numFace,
                         int numDOF);

  // Same as above, but mark the DOF as ESSENTIAL.
  // Multiple degrees of freedom defined on the same entity (vertex/element/edge) cannot
  // be part UNKNOWN and part ESSENTIAL: the whole entity (and thus all of its DOF) is marked as
  // ESSENTIAL.
  void setEssentialVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numVertex,
                             int numDOF = -1);
  void setEssentialElementDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                              int numDOF = -1);
  void setEssentialEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numEdge,
                           int numDOF = -1);
  void setEssentialFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numFace,
                           int numDOF = -1);

  // Get the global number of the specified DOF
  int getVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numVertex,
                   int numDOF = 0);
  int getDOFNumberAtVertex(int iVertex) { return _numberingVertices[iVertex]; };
  int getDOFNumberAtVertex(int iVertex, int numDOF)
  {
    return _numberingVertices[_maxDOFperVertex * iVertex + numDOF];
  };
  int  getElementDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numDOF);
  void setElementDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numDOF, int dofNumber);
  int getEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numEdge, int numDOF);
  int getFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numFace, int numDOF);

  int getDOFCodeAtVertex(int iVertex) { return _codeDOFVertices[iVertex]; };

  bool hasVertexDOFs() { return _maxDOFperVertex > 0; };
  bool hasEdgeDOFs() { return _maxDOFperEdge > 0; };
  bool hasFaceDOFs() { return _maxDOFperFace > 0; };
  bool hasElementDOFs() { return _maxDOFperElem > 0; };

  // Export the vertex DOFs and their coordinates in file.
  // Use the same function in feMetaNumber instead.
  void exportNumberingVertices(feMesh *mesh, FILE *file);

protected:
  void allocateStructures();
  // Set the DOF codes
  void prepareNumbering();
  // Assign the global number to UNKNOWN DOF starting from globalNum.
  // Return the new maximum global number.
  int numberUnknowns(int globalNum);
  // Assign the global number to ESSENTIAL DOF starting from globalNum.
  // Return the new maximum global number.
  // Must be called AFTER numberUnknowns, since ESSENTIAL DOF overwrite
  // the UNKNOWN ones.
  int numberEssential(int globalNum);
  // Fill the allEssential and allUnknown vectors
  void compactFieldDOF();

public:
  void printNumberingVertices(std::ostream &os = std::cout);
  void printNumberingElements(std::ostream &os = std::cout);
  void printNumberingEdges(std::ostream &os = std::cout);
  void printCodeVertices(std::ostream &os = std::cout);
  void printCodeElements(std::ostream &os = std::cout);
  void printCodeEdges(std::ostream &os = std::cout);
};

//
// Handles the numbering of the degrees of freedom (DOF) for all fields.
//
class feMetaNumber
{
protected:
  // Number of unknown DOF for all fields
  int _nInc;
  // Total number of DOF for all fields
  int _nDofs;
  // Number of fields to number
  int _nFields;
  // Name of the fields
  std::vector<std::string> _fieldIDs;

  // Field numberings: <fieldID : feNumber>
  std::map<std::string, feNumber *> _numberings;

public:
  // Create a global numbering for the degrees of freedom of the problem.
  // spaces is the vector of all the FE spaces defined on the mesh, and
  // essentialSpaces is the vector of FE spaces on which essential
  // boundary conditions are applied.
  feMetaNumber(feMesh *mesh,
               const std::vector<feSpace*> &spaces,
               const std::vector<feSpace*> &essentialSpaces);
  feMetaNumber(){};
  ~feMetaNumber();

  int getNbUnknowns() { return _nInc; }
  int getNbDOFs() { return _nDofs; }
  int getNbDOFs(std::string fieldID) { return _numberings[fieldID]->getNbDOFs(); }
  int getNbFields() { return _nFields; }
  std::string getFieldID(int iField) { return _fieldIDs[iField]; }

  // Return the numbering of target field
  feNumber *getNumbering(std::string fieldID) { return _numberings[fieldID]; }
  feNumber *getNumbering(int fieldTag) { return _numberings[_fieldIDs[fieldTag]]; }

  // Return a reference to the array of all UNKNOWN DOF (vert/elem/edge mixed)
  // for target field
  std::set<int> &getUnknownDOF(std::string fieldID)
  {
    return _numberings[fieldID]->getUnknownDOF();
  };
  std::set<int> &getUnknownDOF(int fieldTag)
  {
    return _numberings[_fieldIDs[fieldTag]]->getUnknownDOF();
  };

  // Return a reference to the array of all ESSENTIAL DOF (vert/elem/edge mixed)
  // for target field
  std::set<int> &getEssentialDOF(std::string fieldID)
  {
    return _numberings[fieldID]->getEssentialDOF();
  };
  std::set<int> &getEssentialDOF(int fieldTag)
  {
    return _numberings[_fieldIDs[fieldTag]]->getEssentialDOF();
  };

  // Export the vertex DOFs and their coordinates in file.
  feStatus exportNumberingVertices(feMesh *mesh, std::string fileName);

  void printFields(std::ostream &os = std::cout);
  void printNumberings(std::ostream &os = std::cout);
  void printCodes(std::ostream &os = std::cout);
};

#endif