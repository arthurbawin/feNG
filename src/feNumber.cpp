#include "feNumber.h"

extern int FE_VERBOSE;
static bool COUNT_ONLY = true;

feNumber::feNumber(feMesh *mesh) : _nNod(mesh->getNumVertices()), _nEdg(mesh->getNumEdges())
{
  _nElm = 0;
  for(auto *cnc : mesh->getCncGeo()){
    _nElm += cnc->getNumElements();
  }

  // A global elem to edge map : if the global (line) elem is also an edge.
  // Edges are numbered starting from 1, so a 0 value means the element is not an edge.
  _elemToEdge.resize(_nElm, 0);
  for(auto const &cnc : mesh->getCncGeo()) {
    for(int iElm = 0; iElm < cnc->getNumElements(); ++iElm) {
      if(cnc->getNumEdgesPerElem() == 1) {
        int iElmGlobal = cnc->getElementConnectivity(iElm);
        _elemToEdge[iElmGlobal] = cnc->getEdgeConnectivity(iElm, 0);
      }
    }
  }

  _nDOFVertices.resize(_nNod);
  _nDOFElements.resize(_nElm);
  _nDOFEdges.resize(_nEdg);
}

void feNumber::setUnknownVertexDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Assigning %d DOF_UNKNOWN at global vertex %d defined on "
      "connectivity %s: local elem %d at vertex %d", 
      numDOF, vert, cncGeoID.data(), numElem, numVertex);
  #endif
  _nDOFVertices[vert] = numDOF;
  if(!COUNT_ONLY){
    for(int i = 0; i < numDOF; ++i)
      _codeDOFVertices[_maxDOFperVertex * vert + i] = DOF_UNKNOWN;
  }
}

void feNumber::setUnknownElementDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Assigning %d DOF_UNKNOWN at global element %d defined on "
      "connectivity %s: local elem %d", 
      numDOF, elem, cncGeoID.data(), numElem);
  #endif
  _nDOFElements[elem] = numDOF;
  if(!COUNT_ONLY){
    for(int i = 0; i < numDOF; ++i)
      _codeDOFElements[_maxDOFperElem * elem + i] = DOF_UNKNOWN;
  }
}

void feNumber::setUnknownEdgeDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                          int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Assigning %d DOF_UNKNOWN at global edge %d defined on "
      "connectivity %s: local elem %d at edge %d", 
      numDOF, edge, cncGeoID.data(), numElem, numEdge);
  #endif
  _nDOFEdges[edge] = numDOF;
  if(!COUNT_ONLY){
    for(int i = 0; i < numDOF; ++i)
      _codeDOFEdges[_maxDOFperEdge * edge + i] = DOF_UNKNOWN;
  }
}

void feNumber::setEssentialVertexDOF(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                        int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Setting DOF_ESSENTIAL at global vertex %d defined on "
      "connectivity %s: local elem %d at vertex %d", 
      vert, cncGeoID.data(), numElem, numVertex);
  #endif
  if(numDOF == -1){
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFVertices[vert]; ++i)
      _codeDOFVertices[_maxDOFperVertex * vert + i] = DOF_ESSENTIAL;
  } else{
    _codeDOFVertices[_maxDOFperVertex * vert + numDOF] = DOF_ESSENTIAL;
  }
}

void feNumber::setEssentialElementDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Setting DOF_ESSENTIAL at global element %d defined on "
      "connectivity %s: local elem %d", 
      elem, cncGeoID.data(), numElem);
  #endif
  if(numDOF == -1){
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFElements[elem]; ++i)
      _codeDOFElements[_maxDOFperElem * elem + i] = DOF_ESSENTIAL;
  } else{
    _codeDOFElements[_maxDOFperElem * elem + numDOF] = DOF_ESSENTIAL;
  }
}

void feNumber::setEssentialEdgeDOF(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                      int numEdge, int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  #if defined(FENG_DEBUG)
    feInfoCond(FE_VERBOSE > 1, "Setting DOF_ESSENTIAL at global edge %d defined on "
      "connectivity %s: local elem %d at edge %d", 
      edge, cncGeoID.data(), numElem, numEdge);
  #endif
  if(numDOF == -1){
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFEdges[edge]; ++i)
      _codeDOFEdges[_maxDOFperEdge * edge + i] = DOF_ESSENTIAL;
  } else{
    _codeDOFEdges[_maxDOFperEdge * edge + numDOF] = DOF_ESSENTIAL;
  }
}

int feNumber::getVertexDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // return _numberingVertices[vert];
  return _numberingVertices[_maxDOFperVertex * vert + numDOF];
}

int feNumber::getElementDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  return _numberingElements[_maxDOFperElem * elem + numDOF];
}

int feNumber::getEdgeDOF(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                         int numDOF)
{
  // In connecEdges, edges are numbered starting from 1 and can be negative
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  return _numberingEdges[_maxDOFperEdge * edge + numDOF];
}

void feNumber::allocateStructures()
{
  // Initialize the numbering arrays
  _maxDOFperVertex = *std::max_element(_nDOFVertices.begin(), _nDOFVertices.end());
  _maxDOFperElem = *std::max_element(_nDOFElements.begin(), _nDOFElements.end());
  _maxDOFperEdge = 0;
  if(_nEdg > 0) {
    _maxDOFperEdge = *std::max_element(_nDOFEdges.begin(), _nDOFEdges.end());
  }
  _codeDOFVertices.resize(_nNod * _maxDOFperVertex, DOF_NOT_ASSIGNED);
  _codeDOFElements.resize(_nElm * _maxDOFperElem, DOF_NOT_ASSIGNED);
  _codeDOFEdges.resize(_nEdg * _maxDOFperEdge, DOF_NOT_ASSIGNED);
  _numberingVertices.resize(_nNod * _maxDOFperVertex, DOF_NOT_ASSIGNED);
  _numberingElements.resize(_nElm * _maxDOFperElem, DOF_NOT_ASSIGNED);
  _numberingEdges.resize(_nEdg * _maxDOFperEdge, DOF_NOT_ASSIGNED);
}

void feNumber::prepareNumbering()
{
  // Copy the code of the DOF into the numbering array. If multiple DOFs are associated
  // to an entity, they all have the same code, i.e. a vertex, element or edge cannot be
  // part unknown and part essential.
  for(int iVertex = 0; iVertex < _nNod; ++iVertex) {
    for(int iDOF = 0; iDOF < _nDOFVertices[iVertex]; ++iDOF){
      // if(_nDOFVertices[_maxDOFperVertex * iVertex + iDOF] >= 1){
        _numberingVertices[_maxDOFperVertex * iVertex + iDOF] = _codeDOFVertices[_maxDOFperVertex * iVertex + iDOF];
      // }
    }
  }

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
      _numberingElements[_maxDOFperElem * iElm + iDOF] = _codeDOFElements[_maxDOFperElem * iElm + iDOF];
  }

  for(int iEdg = 0; iEdg < _nEdg; ++iEdg) {
    for(int iDOF = 0; iDOF < _nDOFEdges[iEdg]; ++iDOF) {
      _numberingEdges[_maxDOFperEdge * iEdg + iDOF] = _codeDOFEdges[_maxDOFperEdge * iEdg + iDOF];
    }
  }
}

int feNumber::numberUnknowns(int globalNum)
{
  _nInc = 0;

  // Number the unknown vertices DOF
  for(int i = 0; i < _nNod; ++i) {
    for(int j = 0; j < _nDOFVertices[i]; ++j) {
      if(_numberingVertices[_maxDOFperVertex * i + j] == DOF_UNKNOWN) {
        ++_nInc;
        _numberingVertices[_maxDOFperVertex * i + j] = globalNum++;
      }
    }
  }

  // Number the unknown elements DOF
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nDOFElements[i]; ++j) {
      if(_numberingElements[_maxDOFperElem * i + j] == DOF_UNKNOWN) {
        ++_nInc;
        _numberingElements[_maxDOFperElem * i + j] = globalNum;
        // Update _numberingEdges : element DOFs may have already been numbered by boundary elements
        if(_elemToEdge[i] != 0) {
          int iEdge = fabs(_elemToEdge[i]) - 1;
          // It should be OK to use j because we should have _nDOFEdges[iEdge] == _nDOFElements[i]
          // if the element is also an edge
          _numberingEdges[_maxDOFperEdge * iEdge + j] = globalNum++;

        } else {
          globalNum++;
        }
      }
    }
  }

  // Number the unknown edges DOF
  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _nDOFEdges[i]; ++j) {
      if(_numberingEdges[_maxDOFperEdge * i + j] == DOF_UNKNOWN) {
        ++_nInc;
        _numberingEdges[_maxDOFperEdge * i + j] = globalNum++;
      }
    }
  }

  return globalNum;
}

int feNumber::numberEssential(int globalNum)
{
  _nDofs = _nInc;

  // Number the essential vertices DOF
  for(int i = 0; i < _nNod; ++i) {
    for(int j = 0; j < _nDOFVertices[i]; ++j) {
      if(_numberingVertices[_maxDOFperVertex * i + j] == DOF_ESSENTIAL) {
        ++_nDofs;
        _numberingVertices[_maxDOFperVertex * i + j] = globalNum++;
      }
    }
  }

  // Number the essential elements DOF
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nDOFElements[i]; ++j) {
      if(_numberingElements[_maxDOFperElem * i + j] == DOF_ESSENTIAL) {
        ++_nDofs;
        _numberingElements[_maxDOFperElem * i + j] = globalNum;
        // Update _numberingEdges : element DOFs may have already been numbered by boundary elements
        if(_elemToEdge[i] != 0) {
          int iEdge = fabs(_elemToEdge[i]) - 1;
          // It should be OK to use j because we should have _nDOFEdges[iEdge] == _nDOFElements[i]
          // if the element is also an edge
          _numberingEdges[_maxDOFperEdge * iEdge + j] = globalNum++;
        } else {
          globalNum++;
        }
      }
    }
  }

  // Number the essential edges DOF
  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _nDOFEdges[i]; ++j) {
      if(_numberingEdges[_maxDOFperEdge * i + j] == DOF_ESSENTIAL) {
        ++_nDofs;
        _numberingEdges[_maxDOFperEdge * i + j] = globalNum++;
      }
    }
  }

  return globalNum;
}

// Is used in feMetaNumber::exportNumberingVertices.
// Writes the DOFs numbering along with the physical coordinates of the DOFs for this field.
// Writes a file with format : (1 if essential 0 otherwise) #dof x_dof y_dof z_dof
void feNumber::exportNumberingVertices(feMesh *mesh, FILE *file)
{
  std::vector<Vertex> &vertices = mesh->getVertices();
  int dofNumber;
  Vertex v;
  for(size_t i = 0; i < vertices.size(); ++i){
    dofNumber = _numberingVertices[i];
    v = vertices[i];
    fprintf(file, "%d %d %+-1.16e %+-1.16e %+-1.16e\n", _codeDOFVertices[i] == DOF_ESSENTIAL, dofNumber, v.x(), v.y(), v.z());
  }
}

void feNumber::compactFieldDOF()
{
  _allEssentialDOF.clear();
  _allUnknownDOF.clear();
  for(size_t i = 0; i < _nNod; ++i) {
    for(int j = 0; j < _maxDOFperVertex; ++j) {
      if(_codeDOFVertices[_maxDOFperVertex * i + j] == DOF_ESSENTIAL){
        _allEssentialDOF.insert(_numberingVertices[_maxDOFperVertex * i + j]);
      } else if(_codeDOFVertices[_maxDOFperVertex * i + j] == DOF_UNKNOWN){
        _allUnknownDOF.insert(_numberingVertices[_maxDOFperVertex * i + j]);
      } else{
        // "Ghost" degrees of freedom are not kept
      }
    }
  }

  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _maxDOFperElem; ++j) {
      if(_codeDOFElements[_maxDOFperElem * i + j] == DOF_ESSENTIAL){
        _allEssentialDOF.insert(_numberingElements[_maxDOFperElem * i + j]);
      } else if(_codeDOFElements[_maxDOFperElem * i + j] == DOF_UNKNOWN){
        _allUnknownDOF.insert(_numberingElements[_maxDOFperElem * i + j]);
      }
    }
  }

  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[_maxDOFperEdge * i + j] == DOF_ESSENTIAL){
        _allEssentialDOF.insert(_numberingEdges[_maxDOFperEdge * i + j]);
      } else if(_codeDOFEdges[_maxDOFperEdge * i + j] == DOF_UNKNOWN){
        _allUnknownDOF.insert(_numberingEdges[_maxDOFperEdge * i + j]);
      }
    }
  }
}

void feNumber::printNumberingVertices()
{
  feInfo("\t Number of vertices = %d", _nNod);
  feInfo("\t Max number of DOF per vertex = %d", _maxDOFperVertex);
  for(int i = 0; i < _nNod; ++i) {
    for(int j = 0; j < _maxDOFperVertex; ++j) {
      if(_codeDOFVertices[i] == DOF_ESSENTIAL)
        feInfo("\t DOF: %d \t %s", _numberingVertices[_maxDOFperVertex * i + j], "DOF_ESSENTIAL");
      else if(_codeDOFVertices[i] == DOF_UNKNOWN)
        feInfo("\t DOF: %d \t %s", _numberingVertices[_maxDOFperVertex * i + j], "DOF_UNKNOWN");
      else if(_codeDOFVertices[i] == DOF_NOT_ASSIGNED)
        feInfo("\t DOF: %d \t %s", _numberingVertices[_maxDOFperVertex * i + j], "DOF_NOT_ASSIGNED");
      else{
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for vertex DOF %d - %d :/", _numberingVertices[_maxDOFperVertex * i + j], i, j);
        exit(-1);
      }
    }
  }
}

void feNumber::printNumberingElements()
{
  feInfo("\t Number of elements = %d", _nElm);
  feInfo("\t Max number of DOF per element = %d", _maxDOFperElem);
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _maxDOFperElem; ++j) {
      if(_codeDOFElements[i] == DOF_ESSENTIAL)
        feInfo("\t elm %d - dof %d : %d \t %s", i, j, _numberingElements[_maxDOFperElem * i + j], "DOF_ESSENTIAL");
      else if(_codeDOFElements[i] == DOF_UNKNOWN)
        feInfo("\t elm %d - dof %d : %d \t %s", i, j, _numberingElements[_maxDOFperElem * i + j], "DOF_UNKNOWN");
      else if(_codeDOFElements[i] == DOF_NOT_ASSIGNED)
        feInfo("\t elm %d - dof %d : %d \t %s", i, j, _numberingElements[_maxDOFperElem * i + j], "DOF_NOT_ASSIGNED");
      else{
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for element DOF %d - %d :/", _numberingElements[_maxDOFperElem * i + j], i, j);
        exit(-1);
      }
    }
  }
}

void feNumber::printNumberingEdges()
{
  feInfo("\t Number of edges = %d", _nEdg);
  feInfo("\t Max number of DOF per edge = %d", _maxDOFperEdge);
  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[i] == DOF_ESSENTIAL)
        feInfo("\t %d \t %s", _numberingEdges[_maxDOFperEdge * i + j], "DOF_ESSENTIAL");
      else if(_codeDOFEdges[i] == DOF_UNKNOWN)
        feInfo("\t %d \t %s", _numberingEdges[_maxDOFperEdge * i + j], "DOF_UNKNOWN");
      else if(_codeDOFEdges[i] == DOF_NOT_ASSIGNED)
        feInfo("\t %d \t %s", _numberingEdges[_maxDOFperEdge * i + j], "DOF_NOT_ASSIGNED");
      else{
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for edge DOF %d - %d :/", _numberingEdges[_maxDOFperEdge * i + j], i, j);
        exit(-1);
      }
    }
  }
}

void feNumber::printCodeVertices()
{
  for(auto val : _codeDOFVertices) std::cout << val << std::endl;
}

void feNumber::printCodeElements()
{
  for(auto val : _codeDOFElements) std::cout << val << std::endl;
}

void feNumber::printCodeEdges()
{
  for(auto val : _codeDOFEdges) std::cout << val << std::endl;
}

feMetaNumber::feMetaNumber(feMesh *mesh, const std::vector<feSpace *> &spaces,
                           const std::vector<feSpace *> &essentialSpaces)
{
  for(auto *fS : spaces){
    if(fS == nullptr){
      feErrorMsg(FE_STATUS_ERROR, "Null pointer in vector of FE spaces, maybe you forgot to initialize it.");
      exit(-1);
    }
  }
  for(auto *fS : essentialSpaces){
    if(fS == nullptr){
      feErrorMsg(FE_STATUS_ERROR, "Null pointer in vector of FE essentialSpaces, maybe you forgot to initialize it.");
      exit(-1);
    }
  }

  for(feSpace *fS : spaces) {
    if(std::find(_fieldIDs.begin(), _fieldIDs.end(), fS->getFieldID()) == _fieldIDs.end())
      _fieldIDs.push_back(fS->getFieldID());
  }
  _nFields = _fieldIDs.size();
  
  for(feSpace *fSBC : essentialSpaces) {
    bool err = true;
    for(feSpace *fS : spaces)
      if(fSBC->getFieldID() == fS->getFieldID() && fSBC->getCncGeoID() == fS->getCncGeoID())
        err = false;
    if(err){
      feErrorMsg(FE_STATUS_ERROR, "Space for field \"%s\" and connectivity \"%s\" was marked as\n"
        " essential but is missing from the vector of FE spaces.\n"
        " All essential spaces must be specified in the vector of FE spaces.",
        fSBC->getFieldID().data(),
        fSBC->getCncGeoID().data());
      exit(-1);
    }
  }

  // Create one numbering for each field
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]] = new feNumber(mesh);
  }

  for(feSpace *fS : spaces) {
    fS->setNumberingPtr(_numberings[fS->getFieldID()]);
  }

  // Count DOFs: first pass with initializeUnknowns
  // (it's a little hack to not have to modify all FE spaces
  // and add a count function)
  COUNT_ONLY = true;

  for(feSpace *fS : spaces) {
    fS->initializeNumberingUnknowns();
  }

  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->allocateStructures();
  }

  COUNT_ONLY = false;

  // Now mark the DOFs as unknown or essential
  for(feSpace *fS : spaces) {
    fS->initializeNumberingUnknowns();
  }
  for(feSpace *fSBC : essentialSpaces) {
    fSBC->initializeNumberingEssential();
  }

  // Copy the codes into the numberingArray
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->prepareNumbering();
  }
  
  // Number the DOFs
  int globalNum = 0;
  for(int i = 0; i < _nFields; ++i) {
    globalNum = _numberings[_fieldIDs[i]]->numberUnknowns(globalNum);
  }
  // printNumberings();
  _nInc = globalNum;
  for(int i = 0; i < _nFields; ++i) {
    globalNum = _numberings[_fieldIDs[i]]->numberEssential(globalNum);
  }
  _nDofs = globalNum;

  // Summarize all the degrees of freedom assigned for each field in the _allFieldDOF vector
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->compactFieldDOF();
  }

  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "DEGREES OF FREEDOM:");
  for(int i = 0; i < _nFields; ++i) {
    feInfoCond(FE_VERBOSE > 0, "\t\tInfo for field \"%s\":",
      _fieldIDs[i].data());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of DOF: %d",
      _numberings[_fieldIDs[i]]->getNbDOFs());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of unknown DOF: %d",
      _numberings[_fieldIDs[i]]->getUnknownDOF().size());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of essential DOF: %d",
      _numberings[_fieldIDs[i]]->getEssentialDOF().size());
  }
}

feMetaNumber::~feMetaNumber()
{
  // Delete all numberings
  for(auto it = _numberings.begin(); it != _numberings.end(); )
  {
    delete it->second;
    it = _numberings.erase(it);
  }
}

feStatus feMetaNumber::exportNumberingVertices(feMesh *mesh, std::string fileName)
{
  FILE *file = fopen(fileName.data(), "w");
  if(file == nullptr){
    return feErrorMsg(FE_STATUS_WRITE_ERROR, "Could not export field numbering to file %s", fileName.data());
  }

  for(auto pair : _numberings){
    pair.second->exportNumberingVertices(mesh, file);
  }

  fclose(file);
  return FE_STATUS_OK;
}

void feMetaNumber::printFields()
{
  for(auto s : _fieldIDs) std::cout << s << std::endl;
}

void feMetaNumber::printNumberings()
{
  for(int i = 0; i < _nFields; ++i) {
    feInfo("Field \"%s\" - Numbering of vertex degrees of freedom:", _fieldIDs[i].data());
    _numberings[_fieldIDs[i]]->printNumberingVertices();
    feInfo("Field \"%s\" - Numbering of element degrees of freedom:", _fieldIDs[i].data());
    _numberings[_fieldIDs[i]]->printNumberingElements();
    feInfo("Field \"%s\" - Numbering of edge degrees of freedom:", _fieldIDs[i].data());
    _numberings[_fieldIDs[i]]->printNumberingEdges();
  }
}

void feMetaNumber::printCodes()
{ // Attention : ordre dépend du nom des variables (ordre du map.first)
  for(auto const &secondIsNumber : _numberings) {
    std::cout << "Field " << secondIsNumber.first << " - vertices :" << std::endl;
    secondIsNumber.second->printCodeVertices();
    std::cout << "Field " << secondIsNumber.first << " - elements :" << std::endl;
    secondIsNumber.second->printCodeElements();
    std::cout << "Field " << secondIsNumber.first << " - edges :" << std::endl;
    secondIsNumber.second->printCodeEdges();
  }
}