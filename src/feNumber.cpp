#include "feNumber.h"

feNumber::feNumber(feMesh *mesh) : _nNod(mesh->getNbNodes()), _nEdg(mesh->getNbEdges())
{
  _nElm = 0;
  for(size_t i = 0; i < mesh->getCncGeo().size(); ++i) _nElm += mesh->getCncGeo()[i]->getNbElm();

  // A global elem to edge map : if the global (line) elem is also an edge.
  // Edges are numbered starting from 1, so a 0 value means the element is not an edge.
  _elemToEdge.resize(_nElm, 0);
  for(auto const &cnc : mesh->getCncGeo()) {
    for(int iElm = 0; iElm < cnc->getNbElm(); ++iElm) {
      if(cnc->getNbEdgePerElem() == 1) {
        int iElmGlobal = cnc->getElementConnectivity(iElm);
        _elemToEdge[iElmGlobal] = cnc->getEdgeConnectivity(iElm, 0);
      }
    }
  }

  _nDOFVertices.resize(_nNod);
  _nDOFElements.resize(_nElm);
  _nDOFEdges.resize(_nEdg);
  _codeDOFVertices.resize(_nNod, DOF_NOT_ASSIGNED);
  _codeDOFElements.resize(_nElm, DOF_NOT_ASSIGNED);
  _codeDOFEdges.resize(_nEdg, DOF_NOT_ASSIGNED);

  // TODO : Ajouter une b-rep et verifier la compatibilite des feSpace frontieres
}

void feNumber::defDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // printf("Assigning DOF_UNKNOWN and 1 dof at vertex %d on cnc %s on elem %d at vertex %d\n", vert,
  // cncGeoID.c_str(), numElem, numVertex);
  _nDOFVertices[vert] = numDOF;
  _codeDOFVertices[vert] = DOF_UNKNOWN;
}

void feNumber::defDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  // printf("Assigning DOF_UNKNOWN and %d dof(s) at elem %d on cnc %s on elem %d\n", numDOF, elem,
  // cncGeoID.c_str(), numElem);
  _nDOFElements[elem] = numDOF;
  _codeDOFElements[elem] = DOF_UNKNOWN;
}

void feNumber::defDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                          int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1; // fabs ?
  // printf("Assigning DOF_UNKNOWN and %d dof(s) at edge %d which is edge number %d of elem %d on cnc %s\n",
  // numDOF, edge, numEdge, numElem, cncGeoID.c_str());
  _nDOFEdges[edge] = numDOF;
  _codeDOFEdges[edge] = DOF_UNKNOWN;
}

void feNumber::defDDLSommet_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                        int numVertex)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // printf("Setting global vertex %d as DOF_ESSENTIAL - on cnc %s on elem %d at vertex %d\n", vert,
  // cncGeoID.c_str(), numElem, numVertex);
  _codeDOFVertices[vert] = DOF_ESSENTIAL;
}

void feNumber::defDDLElement_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  // printf("Setting global element %d as DOF_ESSENTIAL - on cnc %s on elem %d\n", elem, cncGeoID.c_str(),
  // numElem);
  _codeDOFElements[elem] = DOF_ESSENTIAL;
}

void feNumber::defDDLEdge_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                      int numEdge)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  // printf("Setting global edge %d as DOF_ESSENTIAL which is edge number %d of elem %d on cnc %s\n", edge,
  // numEdge, numElem, cncGeoID.c_str());
  _codeDOFEdges[edge] = DOF_ESSENTIAL;
}

int feNumber::getDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // return _numberingVertices[vert];
  return _numberingVertices[_maxDOFperVertex * vert + numDOF];
}

int feNumber::getDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  return _numberingElements[_maxDOFperElem * elem + numDOF];
}

int feNumber::getDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                         int numDOF)
{
  // In connecEdges, edges are numbered starting from 1 and can be negative
  int edge = abs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  return _numberingEdges[_maxDOFperEdge * edge + numDOF];
}

void feNumber::prepareNumbering()
{
  _maxDOFperVertex = *std::max_element(_nDOFVertices.begin(), _nDOFVertices.end());
  _maxDOFperElem = *std::max_element(_nDOFElements.begin(), _nDOFElements.end());
  _maxDOFperEdge = 0;
  if(_nEdg > 0) {
    _maxDOFperEdge = *std::max_element(_nDOFEdges.begin(), _nDOFEdges.end());
  }

  // Initialize the numbering arrays
  _numberingVertices.resize(_nNod * _maxDOFperVertex, DOF_NOT_ASSIGNED);
  _numberingElements.resize(_nElm * _maxDOFperElem, DOF_NOT_ASSIGNED);
  _numberingEdges.resize(_nEdg * _maxDOFperEdge, DOF_NOT_ASSIGNED);

  // Copy the code of the DOF into the numbering array. If multiple DOFs are associated
  // to an entity, they all have the same code, i.e. a vertex, element or edge cannot be
  // part unknown and part essential.
  for(int iVertex = 0; iVertex < _nNod; ++iVertex) {
    for(int iDOF = 0; iDOF < _nDOFVertices[iVertex]; ++iDOF){
      // if(_nDOFVertices[_maxDOFperVertex * iVertex + iDOF] >= 1){
        _numberingVertices[_maxDOFperVertex * iVertex + iDOF] = _codeDOFVertices[iVertex];
      // }
    }
  }

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
      _numberingElements[_maxDOFperElem * iElm + iDOF] = _codeDOFElements[iElm];
  }

  for(int iEdg = 0; iEdg < _nEdg; ++iEdg) {
    for(int iDOF = 0; iDOF < _nDOFEdges[iEdg]; ++iDOF) {
      _numberingEdges[_maxDOFperEdge * iEdg + iDOF] = _codeDOFEdges[iEdg];
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
  for(size_t i = 0; i < _numberingVertices.size(); ++i) {
    if(_codeDOFVertices[i] == DOF_ESSENTIAL){
      _allEssentialDOF.push_back(_numberingVertices[i]);
    } else if(_codeDOFVertices[i] == DOF_UNKNOWN){
      _allUnknownDOF.push_back(_numberingVertices[i]);
    } else{
      // "Ghost" degrees of freedom are not kept
    }
  }

  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _maxDOFperElem; ++j) {
      if(_codeDOFElements[i] == DOF_ESSENTIAL){
        _allEssentialDOF.push_back(_numberingElements[_maxDOFperElem * i + j]);
      } else if(_codeDOFElements[i] == DOF_UNKNOWN){
        _allUnknownDOF.push_back(_numberingElements[_maxDOFperElem * i + j]);
      }
    }
  }

  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[i] == DOF_ESSENTIAL){
        _allEssentialDOF.push_back(_numberingEdges[_maxDOFperEdge * i + j]);
      } else if(_codeDOFEdges[i] == DOF_UNKNOWN){
        _allUnknownDOF.push_back(_numberingEdges[_maxDOFperEdge * i + j]);
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
  // AJOUTECHAMP
  for(feSpace *fS : spaces) {
    if(std::find(_fieldIDs.begin(), _fieldIDs.end(), fS->getFieldID()) == _fieldIDs.end())
      _fieldIDs.push_back(fS->getFieldID());
  }
  _nFields = _fieldIDs.size();
  // VRFCHAMP
  for(feSpace *fSBC : essentialSpaces) {
    bool err = true;
    for(feSpace *fS : spaces)
      if(fSBC->getFieldID() == fS->getFieldID() && fSBC->getCncGeoID() == fS->getCncGeoID())
        err = false;
    if(err)
      feWarning("Condition essentielle imposée sur un champ ou une connectivité non présent(e).");
  }
  // Une numerotation pour chaque champ
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]] = new feNumber(mesh);
  }

  // INITIALISE_FENUMER
  for(feSpace *fS : spaces) {
    fS->initializeNumberingUnknowns(_numberings[fS->getFieldID()]);
  }
  // INITIALISE_FENUMER_CLMDOF_ESSENTIAL
  for(feSpace *fSBC : essentialSpaces) {
    fSBC->initializeNumberingEssential(_numberings[fSBC->getFieldID()]);
  }
  // PREPARER_NUMEROTATION
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->prepareNumbering();
  }
  // NUMEROTER
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

}

feMetaNumber::~feMetaNumber()
{
  // Delete all numberings
  for(std::map<std::string, feNumber *>::iterator it = _numberings.begin(); it != _numberings.end(); )
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