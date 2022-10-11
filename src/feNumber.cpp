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
  _codeDOFVertices.resize(_nNod, -3);
  _codeDOFElements.resize(_nElm, -3);
  _codeDOFEdges.resize(_nEdg, -3);

  // TODO : Ajouter une b-rep et verifier la compatibilite des feSpace frontieres
}

void feNumber::defDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // printf("Assigning INC and 1 dof at vertex %d on cnc %s on elem %d at vertex %d\n", vert,
  // cncGeoID.c_str(), numElem, numVertex);
  _nDOFVertices[vert] = 1;
  _codeDOFVertices[vert] = INC;
}

void feNumber::defDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  // printf("Assigning INC and %d dof(s) at elem %d on cnc %s on elem %d\n", numDOF, elem,
  // cncGeoID.c_str(), numElem);
  _nDOFElements[elem] = numDOF;
  _codeDOFElements[elem] = INC;
}

void feNumber::defDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                          int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1; // fabs ?
  // printf("Assigning INC and %d dof(s) at edge %d which is edge number %d of elem %d on cnc %s\n",
  // numDOF, edge, numEdge, numElem, cncGeoID.c_str());
  _nDOFEdges[edge] = numDOF;
  _codeDOFEdges[edge] = INC;
}

void feNumber::defDDLSommet_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                        int numVertex)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  // printf("Setting global vertex %d as ESS - on cnc %s on elem %d at vertex %d\n", vert,
  // cncGeoID.c_str(), numElem, numVertex);
  _codeDOFVertices[vert] = ESS;
}

void feNumber::defDDLElement_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  // printf("Setting global element %d as ESS - on cnc %s on elem %d\n", elem, cncGeoID.c_str(),
  // numElem);
  _codeDOFElements[elem] = ESS;
}

void feNumber::defDDLEdge_essentialBC(feMesh *mesh, std::string const &cncGeoID, int numElem,
                                      int numEdge)
{
  int edge = fabs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  // printf("Setting global edge %d as ESS which is edge number %d of elem %d on cnc %s\n", edge,
  // numEdge, numElem, cncGeoID.c_str());
  _codeDOFEdges[edge] = ESS;
}

int feNumber::getDDLSommet(feMesh *mesh, std::string const &cncGeoID, int numElem, int numVertex)
{
  int vert = mesh->getVertex(cncGeoID, numElem, numVertex);
  return _numberingVertices[vert];
}

int feNumber::getDDLElement(feMesh *mesh, std::string const &cncGeoID, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoID, numElem);
  return _numberingElements[_maxDOFperElem * elem + numDOF];
}

int feNumber::getDDLEdge(feMesh *mesh, std::string const &cncGeoID, int numElem, int numEdge,
                         int numDOF)
{
  // In connecEdges, edges are numbering starting from 1 and can be negative
  int edge = abs(mesh->getEdge(cncGeoID, numElem, numEdge)) - 1;
  return _numberingEdges[_maxDOFperEdge * edge + numDOF];
}

void feNumber::prepareNumbering()
{
  _maxDOFperElem = *std::max_element(_nDOFElements.begin(), _nDOFElements.end());
  _maxDOFperEdge = 0;
  if(_nEdg > 0) {
    _maxDOFperEdge = *std::max_element(_nDOFEdges.begin(), _nDOFEdges.end());
  }

  _numberingVertices.resize(_nNod, -3);
  _numberingElements.resize(_nElm * _maxDOFperElem, -3);
  _numberingEdges.resize(_nEdg * _maxDOFperEdge, -3);

  for(int iVertex = 0; iVertex < _nNod; ++iVertex) {
    if(_nDOFVertices[iVertex] == 1) _numberingVertices[iVertex] = _codeDOFVertices[iVertex];
  }
  // printf("Numérotation des sommets : \n");
  // for(int iVertex = 0; iVertex < _nNod; ++iVertex) {
  //   std::cout<<_numberingVertices[iVertex]<<std::endl;
  // }
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
      _numberingElements[_maxDOFperElem * iElm + iDOF] = _codeDOFElements[iElm];
  }
  // printf("Numérotation des elements : \n");
  // for(int iElm = 0; iElm < _nElm; ++iElm) {
  //   for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
  //     std::cout<<_numberingElements[_maxDOFperElem * iElm + iDOF]<<std::endl;
  // }
  for(int iEdg = 0; iEdg < _nEdg; ++iEdg) {
    for(int iDOF = 0; iDOF < _nDOFEdges[iEdg]; ++iDOF) {
      _numberingEdges[_maxDOFperEdge * iEdg + iDOF] = _codeDOFEdges[iEdg];
    }
  }
  // printf("Numérotation des aretes : \n");
  // for(int iEdg = 0; iEdg < _nEdg; ++iEdg) {
  //   for(int iDOF = 0; iDOF < _nDOFEdges[iEdg]; ++iDOF) {
  //     std::cout<<_numberingEdges[_maxDOFperEdge * iEdg + iDOF]<<std::endl;
  //   }
  // }
}

int feNumber::numberUnknowns(int globalNum)
{
  _nInc = 0;
  for(int i = 0; i < _nNod; ++i) {
    if(_numberingVertices[i] == INC) {
      ++_nInc;
      _numberingVertices[i] = globalNum++;
    }
  }
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nDOFElements[i]; ++j) {
      if(_numberingElements[_maxDOFperElem * i + j] == INC) {
        ++_nInc;
        _numberingElements[_maxDOFperElem * i + j] = globalNum;
        // std::cout<<"coucou1"<<std::endl;
        // Update _numberingEdges : element DOFs may have already been numbered by boundary elements
        if(_elemToEdge[i] != 0) {
          // std::cout<<"coucou2"<<std::endl;
          int iEdge = fabs(_elemToEdge[i]) - 1;
          // std::cout<<"coucou3"<<std::endl;
          // It should be OK to use j because we should have _nDOFEdges[iEdge] == _nDOFElements[i]
          // if the element is also an edge
          _numberingEdges[_maxDOFperEdge * iEdge + j] = globalNum++;
          // std::cout<<"coucou4"<<std::endl;

        } else {
          globalNum++;
        }
      }
    }
  }
  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _nDOFEdges[i]; ++j) {
      if(_numberingEdges[_maxDOFperEdge * i + j] == INC) {
        ++_nInc;
        _numberingEdges[_maxDOFperEdge * i + j] = globalNum++;
        // std::cout<<"coucou"<<std::endl;
      }
    }
  }
  return globalNum;
}

int feNumber::numberEssential(int globalNum)
{
  _nDofs = _nInc;
  for(int i = 0; i < _nNod; ++i) {
    if(_numberingVertices[i] == ESS) {
      ++_nDofs;
      _numberingVertices[i] = globalNum++;
    }
  }
  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _nDOFElements[i]; ++j) {
      if(_numberingElements[_maxDOFperElem * i + j] == ESS) {
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
  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _nDOFEdges[i]; ++j) {
      if(_numberingEdges[_maxDOFperEdge * i + j] == ESS) {
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
  double x, y, z;
  Vertex v;
  for(size_t i = 0; i < vertices.size(); ++i){
    dofNumber = _numberingVertices[i];
    v = vertices[i];
    fprintf(file, "%d %d %+-1.16e %+-1.16e %+-1.16e\n", _codeDOFVertices[i] == ESS, dofNumber, v.x(), v.y(), v.z());
  }
}

void feNumber::compactFieldDOF()
{
  _allEssentialDOF.clear();
  _allUnknownDOF.clear();
  for(size_t i = 0; i < _numberingVertices.size(); ++i) {
    if(_codeDOFVertices[i] == ESS){
      _allEssentialDOF.push_back(_numberingVertices[i]);
    } else if(_codeDOFVertices[i] == INC){
      _allUnknownDOF.push_back(_numberingVertices[i]);
    } else{
      // "Ghost" degrees of freedom are not kept
    }
  }

  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _maxDOFperElem; ++j) {
      if(_codeDOFElements[i] == ESS){
        _allEssentialDOF.push_back(_numberingElements[_maxDOFperElem * i + j]);
      } else if(_codeDOFElements[i] == INC){
        _allUnknownDOF.push_back(_numberingElements[_maxDOFperElem * i + j]);
      }
    }
  }

  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[i] == ESS){
        _allEssentialDOF.push_back(_numberingEdges[_maxDOFperEdge * i + j]);
      } else if(_codeDOFEdges[i] == INC){
        _allUnknownDOF.push_back(_numberingEdges[_maxDOFperEdge * i + j]);
      }
    }
  }
}

void feNumber::printNumberingVertices()
{
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

void feNumber::printNumberingElements()
{
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

void feNumber::printNumberingEdges()
{
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
  // INITIALISE_FENUMER_CLMESS
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
    // for(auto const& secondIsNumber : _numberings){
    std::cout << "Field " << _fieldIDs[i] << " - vertices :" << std::endl;
    _numberings[_fieldIDs[i]]->printNumberingVertices();
    std::cout << "Field " << _fieldIDs[i] << " - elements :" << std::endl;
    _numberings[_fieldIDs[i]]->printNumberingElements();
    std::cout << "Field " << _fieldIDs[i] << " - edges :" << std::endl;
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