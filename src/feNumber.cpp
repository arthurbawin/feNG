#include "feNumber.h"

extern int FE_VERBOSE;
static bool COUNT_ONLY = true;

feNumber::feNumber(feMesh *mesh, const std::string &fieldName) :
  _fieldID(fieldName),
  _nNod(mesh->getNumVertices()),
  _nEdg(mesh->getNumEdges()),
  _nFac(mesh->getNumFaces())
{
  _nElm = 0;
  for(auto *cnc : mesh->getCncGeo()) {
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

  // TODO: same map for triangle facets

  _nDOFVertices.resize(_nNod);
  _nDOFElements.resize(_nElm);
  _nDOFEdges.resize(_nEdg);
  _nDOFFaces.resize(_nFac);
}

void feNumber::setUnknownVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                   int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoTag, numElem, numVertex);
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Assigning %d DOF_UNKNOWN at global vertex %d defined on "
             "connectivity %d: local elem %d at vertex %d",
             numDOF, vert, cncGeoTag, numElem, numVertex);
#endif
  _nDOFVertices[vert] = numDOF;
  if(!COUNT_ONLY) {
    for(int i = 0; i < numDOF; ++i) _codeDOFVertices[_maxDOFperVertex * vert + i] = DOF_UNKNOWN;
  }
}

void feNumber::setUnknownElementDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                    int numDOF)
{
  int elem = mesh->getElement(cncGeoTag, numElem);
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Assigning %d DOF_UNKNOWN at global element %d defined on "
             "connectivity %d: local elem %d",
             numDOF, elem, cncGeoTag, numElem);
#endif
  _nDOFElements[elem] = numDOF;
  if(!COUNT_ONLY) {
    for(int i = 0; i < numDOF; ++i) _codeDOFElements[_maxDOFperElem * elem + i] = DOF_UNKNOWN;
  }
}

void feNumber::setUnknownEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                 int numEdge, int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoTag, numElem, numEdge)) - 1;
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Assigning %d DOF_UNKNOWN at global edge %d defined on "
             "connectivity %d: local elem %d at edge %d",
             numDOF, edge, cncGeoTag, numElem, numEdge);
#endif
  _nDOFEdges[edge] = numDOF;
  if(!COUNT_ONLY) {
    for(int i = 0; i < numDOF; ++i) _codeDOFEdges[_maxDOFperEdge * edge + i] = DOF_UNKNOWN;
  }
}

void feNumber::setUnknownFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                 int numFace, int numDOF)
{
  int face = fabs(mesh->getFace(cncGeoTag, numElem, numFace)) - 1;
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Assigning %d DOF_UNKNOWN at global face %d defined on "
             "connectivity %d: local elem %d at face %d",
             numDOF, face, cncGeoTag, numElem, numFace);
#endif
  _nDOFFaces[face] = numDOF;
  if(!COUNT_ONLY) {
    for(int i = 0; i < numDOF; ++i) _codeDOFFaces[_maxDOFperFace * face + i] = DOF_UNKNOWN;
  }
}

void feNumber::setEssentialVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                     int numVertex, int numDOF)
{
  int vert = mesh->getVertex(cncGeoTag, numElem, numVertex);
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Setting DOF_ESSENTIAL at global vertex %d defined on "
             "connectivity %d: local elem %d at vertex %d",
             vert, cncGeoTag, numElem, numVertex);
#endif
  if(numDOF == -1) {
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFVertices[vert]; ++i)
      _codeDOFVertices[_maxDOFperVertex * vert + i] = DOF_ESSENTIAL;
  } else {
    _codeDOFVertices[_maxDOFperVertex * vert + numDOF] = DOF_ESSENTIAL;
  }
}

void feNumber::setEssentialElementDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                      int numDOF)
{
  int elem = mesh->getElement(cncGeoTag, numElem);
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Setting DOF_ESSENTIAL at global element %d defined on "
             "connectivity %d: local elem %d",
             elem, cncGeoTag, numElem);
#endif
  if(numDOF == -1) {
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFElements[elem]; ++i)
      _codeDOFElements[_maxDOFperElem * elem + i] = DOF_ESSENTIAL;
  } else {
    _codeDOFElements[_maxDOFperElem * elem + numDOF] = DOF_ESSENTIAL;
  }
}

void feNumber::setEssentialEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                   int numEdge, int numDOF)
{
  int edge = fabs(mesh->getEdge(cncGeoTag, numElem, numEdge)) - 1;
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Setting DOF_ESSENTIAL at global edge %d defined on "
             "connectivity %d: local elem %d at edge %d",
             edge, cncGeoTag, numElem, numEdge);
#endif
  if(numDOF == -1) {
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFEdges[edge]; ++i)
      _codeDOFEdges[_maxDOFperEdge * edge + i] = DOF_ESSENTIAL;
  } else {
    _codeDOFEdges[_maxDOFperEdge * edge + numDOF] = DOF_ESSENTIAL;
  }
}

void feNumber::setEssentialFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem,
                                   int numFace, int numDOF)
{
  int face = fabs(mesh->getFace(cncGeoTag, numElem, numFace)) - 1;
#if defined(FENG_DEBUG)
  feInfoCond(FE_VERBOSE > 2,
             "Setting DOF_ESSENTIAL at global face %d defined on "
             "connectivity %d: local elem %d at face %d",
             face, cncGeoTag, numElem, numFace);
#endif
  if(numDOF == -1) {
    // Mark all DOFs on this entity as essential
    for(int i = 0; i < _nDOFFaces[face]; ++i)
      _codeDOFFaces[_maxDOFperFace * face + i] = DOF_ESSENTIAL;
  } else {
    _codeDOFFaces[_maxDOFperFace * face + numDOF] = DOF_ESSENTIAL;
  }
}

int feNumber::getVertexDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numVertex,
                           int numDOF)
{
  int vert = mesh->getVertex(cncGeoTag, numElem, numVertex);
  assert(_numberingVertices[_maxDOFperVertex * vert + numDOF] != DOF_NOT_ASSIGNED);
  return _numberingVertices[_maxDOFperVertex * vert + numDOF];
}

int feNumber::getElementDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numDOF)
{
  int elem = mesh->getElement(cncGeoTag, numElem);
  assert(_numberingElements[_maxDOFperElem * elem + numDOF] != DOF_NOT_ASSIGNED);
  return _numberingElements[_maxDOFperElem * elem + numDOF];
}

void feNumber::setElementDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numDOF, int dofNumber)
{
  int elem = mesh->getElement(cncGeoTag, numElem);
  assert(_numberingElements[_maxDOFperElem * elem + numDOF] != DOF_NOT_ASSIGNED);
  _numberingElements[_maxDOFperElem * elem + numDOF] = dofNumber;
}

int feNumber::getEdgeDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numEdge,
                         int numDOF)
{
  // In connecEdges, edges are numbered starting from 1 and can be negative
  int edge = fabs(mesh->getEdge(cncGeoTag, numElem, numEdge)) - 1;
  assert(_numberingEdges[_maxDOFperEdge * edge + numDOF] != DOF_NOT_ASSIGNED);
  return _numberingEdges[_maxDOFperEdge * edge + numDOF];
}

int feNumber::getFaceDOF(feMesh *mesh, const int cncGeoTag, int numElem, int numFace,
                         int numDOF)
{
  // In _connecTriangles, faces are numbered starting from 1 and can be negative
  int face = fabs(mesh->getFace(cncGeoTag, numElem, numFace)) - 1;
  assert(_numberingFaces[_maxDOFperFace * face + numDOF] != DOF_NOT_ASSIGNED);
  return _numberingFaces[_maxDOFperFace * face + numDOF];
}

void feNumber::applyPeriodicity(feMesh *mesh, const std::vector<feSpace*> &spaces)
{
  const double tol = 1e-8;

  for(const auto &s : spaces)
  {
    if(s->getFieldID() == _fieldID && s->isPeriodic() && s->isPeriodicMaster())
    {
      const feSpace *matchingSpace = s->MatchingPeriodicSpace();
      const feCncGeo *cnc0 = mesh->getCncGeoByTag(s->getCncGeoTag());
      const feCncGeo *cnc1 = mesh->getCncGeoByTag(matchingSpace->getCncGeoTag());

      //
      // Vertex DOFs
      //
      const std::vector<int> &vertices0 = cnc0->GlobalVertexTags();
      const std::vector<int> &vertices1 = cnc1->GlobalVertexTags();
      const std::vector<Vertex> &vertices = mesh->getVertices();

      Vertex offset(s->PeriodicOffset());

      int matchingVertices = 0;
      for(size_t i = 0; i < vertices0.size(); ++i)
      {
        const int iv0 = vertices0[i];
        const Vertex &v0 = vertices[iv0];
        Vertex v0translated = v0 + offset;

        for(size_t j = 0; j < vertices1.size(); ++j)
        {
          const int iv1 = vertices1[j];
          const Vertex &v1 = vertices[iv1];
          if(v1.distance(v0translated) < tol)
          {
            // Vertices match: add all their DOF to map
            matchingVertices++;
            for(int iDOF = 0; iDOF < _nDOFVertices[iv0]; ++iDOF)
            {
              const int masterDOF = _numberingVertices[_maxDOFperVertex * iv0 + iDOF];
              const int slaveDOF  = _numberingVertices[_maxDOFperVertex * iv1 + iDOF];
              _periodicDOF.insert({masterDOF, slaveDOF});
            }
            break;
          }
        }
      }

      //
      // Edge DOFs
      //
      int matchingEdges = 0;
      const std::vector<int> &edges0 = cnc0->GlobalEdgeTags();
      const std::vector<int> &edges1 = cnc1->GlobalEdgeTags();
      const std::map<int, const Edge *> &edges = mesh->getEdgesMap();

      for(size_t i = 0; i < edges0.size(); ++i)
      {
        const int iedge0 = fabs(edges0[i]) - 1;
        const Vertex *v00 = edges.at(edges0[i])->getVertex(0);
        const Vertex *v01 = edges.at(edges0[i])->getVertex(1);
        Vertex v00translated = *v00 + offset;
        Vertex v01translated = *v01 + offset;

        for(size_t j = 0; j < edges1.size(); ++j)
        {
          const int iedge1 = fabs(edges1[j]) - 1;
          const Vertex *v10 = edges.at(edges1[j])->getVertex(0);
          const Vertex *v11 = edges.at(edges1[j])->getVertex(1);

          if((v10->distance(v00translated) < tol && v11->distance(v01translated) < tol) ||
             (v11->distance(v00translated) < tol && v10->distance(v01translated) < tol))
          {
            // Edges match: add all their DOF to map
            matchingEdges++;
            for(int iDOF = 0; iDOF < _nDOFEdges[iedge0]; ++iDOF)
            {
              const int masterDOF = _numberingEdges[_maxDOFperEdge * iedge0 + iDOF];
              const int slaveDOF  = _numberingEdges[_maxDOFperEdge * iedge1 + iDOF];
              _periodicDOF.insert({masterDOF, slaveDOF});
            }
            break;         
          }
        }
      }

      if(matchingVertices == 0 && matchingEdges == 0)
      {
        feErrorMsg(FE_STATUS_ERROR,
          "Could not find any matching vertex or edge when assigning\n"
          " periodic degrees of freedom with tolerance %1.3e. Increase tolerance,\n"
          " or check offset vector for inconsistency with the mesh.\n", tol);
        exit(-1);
      }
    }
  }
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
  _maxDOFperFace = 0;
  if(_nFac > 0) {
    _maxDOFperFace = *std::max_element(_nDOFFaces.begin(), _nDOFFaces.end());
  }
  _codeDOFVertices.resize(_nNod * _maxDOFperVertex, DOF_NOT_ASSIGNED);
  _codeDOFElements.resize(_nElm * _maxDOFperElem, DOF_NOT_ASSIGNED);
  _codeDOFEdges.resize(_nEdg * _maxDOFperEdge, DOF_NOT_ASSIGNED);
  _codeDOFFaces.resize(_nFac * _maxDOFperFace, DOF_NOT_ASSIGNED);
  _numberingVertices.resize(_nNod * _maxDOFperVertex, DOF_NOT_ASSIGNED);
  _numberingElements.resize(_nElm * _maxDOFperElem, DOF_NOT_ASSIGNED);
  _numberingEdges.resize(_nEdg * _maxDOFperEdge, DOF_NOT_ASSIGNED);
  _numberingFaces.resize(_nFac * _maxDOFperFace, DOF_NOT_ASSIGNED);
}

void feNumber::prepareNumbering()
{
  // Copy the code of the DOF into the numbering array. If multiple DOFs are associated
  // to an entity, they all have the same code, i.e. a vertex, element or edge cannot be
  // part unknown and part essential.
  for(int iVertex = 0; iVertex < _nNod; ++iVertex) {
    for(int iDOF = 0; iDOF < _nDOFVertices[iVertex]; ++iDOF) {
      // if(_nDOFVertices[_maxDOFperVertex * iVertex + iDOF] >= 1){
      _numberingVertices[_maxDOFperVertex * iVertex + iDOF] =
        _codeDOFVertices[_maxDOFperVertex * iVertex + iDOF];
      // }
    }
  }

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int iDOF = 0; iDOF < _nDOFElements[iElm]; ++iDOF)
      _numberingElements[_maxDOFperElem * iElm + iDOF] =
        _codeDOFElements[_maxDOFperElem * iElm + iDOF];
  }

  for(int iEdg = 0; iEdg < _nEdg; ++iEdg) {
    for(int iDOF = 0; iDOF < _nDOFEdges[iEdg]; ++iDOF) {
      _numberingEdges[_maxDOFperEdge * iEdg + iDOF] = _codeDOFEdges[_maxDOFperEdge * iEdg + iDOF];
    }
  }

  for(int iFace = 0; iFace < _nFac; ++iFace) {
    for(int iDOF = 0; iDOF < _nDOFFaces[iFace]; ++iDOF) {
      _numberingFaces[_maxDOFperFace * iFace + iDOF] = _codeDOFFaces[_maxDOFperFace * iFace + iDOF];
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
  // TODO : Like for the edges, identify the boundary triangles which are also faces of tetrahedra
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

  // Number the unknown faces DOF
  for(int i = 0; i < _nFac; ++i) {
    for(int j = 0; j < _nDOFFaces[i]; ++j) {
      if(_numberingFaces[_maxDOFperFace * i + j] == DOF_UNKNOWN) {
        ++_nInc;
        _numberingFaces[_maxDOFperFace * i + j] = globalNum++;
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
  // TODO : Modify for faces like above
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

  // Number the essential edges DOF
  for(int i = 0; i < _nFac; ++i) {
    for(int j = 0; j < _nDOFFaces[i]; ++j) {
      if(_numberingFaces[_maxDOFperFace * i + j] == DOF_ESSENTIAL) {
        ++_nDofs;
        _numberingFaces[_maxDOFperFace * i + j] = globalNum++;
      }
    }
  }

  return globalNum;
}

// Is used in feMetaNumber::exportNumberingVertices.
// Writes the DOFs numbering along with the physical coordinates of the DOFs for this field.
// Writes a file with format : (1 if essential 0 otherwise) #dof x_dof y_dof z_dof
void feNumber::exportNumberingVertices(feMesh *mesh, FILE *file, bool posFormat)
{
  const std::vector<Vertex> &vertices = mesh->getVertices();

  for(size_t i = 0; i < vertices.size(); ++i)
  {
    const int dofNumber = _numberingVertices[i];
    const Vertex &v = vertices[i];

    if(posFormat)
    {
      fprintf(file, "SP(%g,%g,%g){%d};\n", v.x(), v.y(), v.z(), dofNumber);
    } else 
    {
      fprintf(file, "%d %d %+-1.16e %+-1.16e %+-1.16e\n", _codeDOFVertices[i] == DOF_ESSENTIAL,
              dofNumber, v.x(), v.y(), v.z());
    }
  }

  const std::vector<const Edge *> &edges = mesh->getEdges();

  for(size_t i = 0; i < edges.size(); ++i)
  {
    const int edge = fabs(edges[i]->getTag()) - 1;
    const int dofNumber = _numberingEdges[edge];
    const Vertex *v0 = edges[i]->getVertex(0);
    const Vertex *v1 = edges[i]->getVertex(1);
    Vertex vMid = (*v0 + *v1) * 0.5;

    if(posFormat)
    {
      fprintf(file, "SP(%g,%g,%g){%d};\n", vMid.x(), vMid.y(), vMid.z(), dofNumber);
    } else 
    {
      fprintf(file, "%d %d %+-1.16e %+-1.16e %+-1.16e\n", _codeDOFEdges[i] == DOF_ESSENTIAL,
              dofNumber, vMid.x(), vMid.y(), vMid.z());
    }
  }
}

void feNumber::compactFieldDOF()
{
  _allEssentialDOF.clear();
  _allUnknownDOF.clear();
  for(int i = 0; i < _nNod; ++i) {
    for(int j = 0; j < _maxDOFperVertex; ++j) {
      if(_codeDOFVertices[_maxDOFperVertex * i + j] == DOF_ESSENTIAL) {
        _allEssentialDOF.insert(_numberingVertices[_maxDOFperVertex * i + j]);
      } else if(_codeDOFVertices[_maxDOFperVertex * i + j] == DOF_UNKNOWN) {
        _allUnknownDOF.insert(_numberingVertices[_maxDOFperVertex * i + j]);
      } else {
        // "Ghost" degrees of freedom are not kept
      }
    }
  }

  for(int i = 0; i < _nElm; ++i) {
    for(int j = 0; j < _maxDOFperElem; ++j) {
      if(_codeDOFElements[_maxDOFperElem * i + j] == DOF_ESSENTIAL) {
        _allEssentialDOF.insert(_numberingElements[_maxDOFperElem * i + j]);
      } else if(_codeDOFElements[_maxDOFperElem * i + j] == DOF_UNKNOWN) {
        _allUnknownDOF.insert(_numberingElements[_maxDOFperElem * i + j]);
      }
    }
  }

  for(int i = 0; i < _nEdg; ++i) {
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[_maxDOFperEdge * i + j] == DOF_ESSENTIAL) {
        _allEssentialDOF.insert(_numberingEdges[_maxDOFperEdge * i + j]);
      } else if(_codeDOFEdges[_maxDOFperEdge * i + j] == DOF_UNKNOWN) {
        _allUnknownDOF.insert(_numberingEdges[_maxDOFperEdge * i + j]);
      }
    }
  }
}

void feNumber::printNumberingVertices(std::ostream &os)
{
  os << "Number of vertices = " << _nNod << std::endl;
  os << "Max number of DOF per vertex = " << _maxDOFperVertex << std::endl;
  for(int i = 0; i < _nNod; ++i) {
    os << "Field " << _fieldID << " - vertex " << i << ": ";
    for(int j = 0; j < _maxDOFperVertex; ++j) {
      int code = _codeDOFVertices[_maxDOFperVertex * i + j];
      if(code == DOF_ESSENTIAL)
        os << _numberingVertices[_maxDOFperVertex * i + j] << " - DOF_ESSENTIAL" << "\t";
      else if(code == DOF_UNKNOWN)
        os << _numberingVertices[_maxDOFperVertex * i + j] << " - DOF_UNKNOWN" << "\t";
      else if(code == DOF_NOT_ASSIGNED)
        os << _numberingVertices[_maxDOFperVertex * i + j] << " - DOF_NOT_ASSIGNED" << "\t";
      else {
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for vertex DOF %d - %d :/",
                   _numberingVertices[_maxDOFperVertex * i + j], i, j);
        exit(-1);
      }
    }
    os << std::endl;
  }
}

void feNumber::printNumberingElements(std::ostream &os)
{
  os << "Number of elements = " << _nElm << std::endl;
  os << "Max number of DOF per element = " << _maxDOFperElem << std::endl;
  for(int i = 0; i < _nElm; ++i) {
    os << "Field " << _fieldID << " - element " << i << ": ";
    for(int j = 0; j < _maxDOFperElem; ++j) {
      int code = _codeDOFElements[_maxDOFperElem * i + j];
      if(code == DOF_ESSENTIAL)
        os << _numberingElements[_maxDOFperElem * i + j] << " - DOF_ESSENTIAL" << "\t";
      else if(code == DOF_UNKNOWN)
        os << _numberingElements[_maxDOFperElem * i + j] << " - DOF_UNKNOWN" << "\t";
      else if(code == DOF_NOT_ASSIGNED)
        os << _numberingElements[_maxDOFperElem * i + j] << " - DOF_NOT_ASSIGNED" << "\t";
      else {
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for element DOF %d - %d :/",
                   _numberingElements[_maxDOFperElem * i + j], i, j);
        exit(-1);
      }
    }
    os << std::endl;
  }
}

void feNumber::printNumberingEdges(std::ostream &os)
{
  os << "Number of edges = " << _nEdg << std::endl;
  os << "Max number of DOF per edge = " << _maxDOFperEdge << std::endl;
  for(int i = 0; i < _nEdg; ++i) {
    os << "Field " << _fieldID << " - edge " << i << ": ";
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      if(_codeDOFEdges[i] == DOF_ESSENTIAL)
        os << _numberingEdges[_maxDOFperEdge * i + j] << " - DOF_ESSENTIAL" << "\t";
      else if(_codeDOFEdges[i] == DOF_UNKNOWN)
        os << _numberingEdges[_maxDOFperEdge * i + j] << " - DOF_UNKNOWN" << "\t";
      else if(_codeDOFEdges[i] == DOF_NOT_ASSIGNED)
        os << _numberingEdges[_maxDOFperEdge * i + j] << " - DOF_NOT_ASSIGNED" << "\t";
      else {
        feErrorMsg(FE_STATUS_ERROR, "Unexpected value %d for edge DOF %d - %d :/",
                   _numberingEdges[_maxDOFperEdge * i + j], i, j);
        exit(-1);
      }
    }
    os << std::endl;
  }
}

void feNumber::printCodeVertices(std::ostream &os)
{
  std::cout << _fieldID << std::endl;
  for(int iVert = 0; iVert < _nNod; ++iVert) {
    os << "Field " << _fieldID << " - vertex " << iVert << ": ";
    for(int j = 0; j < _maxDOFperVertex; ++j) {
      os << _codeDOFVertices[iVert * _maxDOFperVertex + j] << "  ";
    }
    os << std::endl;
  }
  // for(auto val : _codeDOFVertices) os << val << std::endl;
}

void feNumber::printCodeEdges(std::ostream &os)
{
  std::cout << _fieldID << std::endl;
  for(int iEdge = 0; iEdge < _nEdg; ++iEdge) {
    os << "Field " << _fieldID << " - edge " << iEdge << ": ";
    for(int j = 0; j < _maxDOFperEdge; ++j) {
      os << _codeDOFEdges[iEdge * _maxDOFperEdge + j] << "  ";
    }
    os << std::endl;
  }
  // for(auto val : _codeDOFEdges) os << val << std::endl;
}

void feNumber::printCodeElements(std::ostream &os)
{
  std::cout << _fieldID << std::endl;
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    os << "Field " << _fieldID << " - element " << iElm << ": ";
    for(int j = 0; j < _maxDOFperElem; ++j) {
      os << _codeDOFElements[iElm * _maxDOFperElem + j] << "  ";
    }
    os << std::endl;
  }
  // for(auto val : _codeDOFElements) os << val << std::endl;
}

// FIXME: should not exit from constructor
// Create a createNumbering wrapper instead for example.
feMetaNumber::feMetaNumber(feMesh *mesh,
                           const std::vector<feSpace*> &spaces,
                           const std::vector<feSpace*> &essentialSpaces)
{
  for(const auto &fS : spaces) {
    if(fS == nullptr) {
      feErrorMsg(FE_STATUS_ERROR,
                 "Null pointer in vector of FE spaces, maybe you forgot to initialize it.");
      exit(-1);
    }
  }
  for(const auto &fS : essentialSpaces) {
    if(fS == nullptr) {
      feErrorMsg(
        FE_STATUS_ERROR,
        "Null pointer in vector of FE essentialSpaces, maybe you forgot to initialize it.");
      exit(-1);
    }
  }

  for(const auto &fS : spaces) {
    if(std::find(_fieldIDs.begin(), _fieldIDs.end(), fS->getFieldID()) == _fieldIDs.end())
      _fieldIDs.push_back(fS->getFieldID());
  }
  _nFields = _fieldIDs.size();

  for(const auto &fSBC : essentialSpaces) {
    bool err = true;
    for(const auto &fS : spaces)
      if(fSBC->representsSameFieldAs(*fS)) err = false;
    if(err) {
      feErrorMsg(FE_STATUS_ERROR,
                 "Space for field \"%s\" and connectivity \"%s\" was marked as\n"
                 " essential but is missing from the vector of FE spaces.\n"
                 " All essential spaces must be specified in the vector of FE spaces.",
                 fSBC->getFieldID().data(), fSBC->getCncGeoID().data());
      exit(-1);
    }
  }

  // Set _essentialComponent to true for all fully essential spaces.
  for(auto &space : essentialSpaces) {
    space->setEssentialComponent(0, true);
    space->setEssentialComponent(1, true);
    space->setEssentialComponent(2, true);
  }

  // Create one numbering for each field
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]] = new feNumber(mesh, _fieldIDs[i]);
  }

  for(auto &fS : spaces) {
    fS->setNumberingPtr(_numberings[fS->getFieldID()]);
  }

  // Count the DOFs: do a first pass with initializeUnknowns
  // and COUNT_ONLY = true (it's a little hack to avoid
  // modifying all FE spaces and adding a count function)
  COUNT_ONLY = true;

  for(const auto &fS : spaces) {
    fS->initializeNumberingUnknowns();
  }

  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->allocateStructures();
  }

  COUNT_ONLY = false;

  // Now mark the DOFs as unknown or essential
  for(const auto &fS : spaces) {
    fS->initializeNumberingUnknowns();
  }
  for(const auto &fSBC : essentialSpaces) {
    fSBC->initializeNumberingEssential();
  }

  // For discontinuous fields: update relevant unknown DOFs
  // which have been marked as essential by another feSpace.
  // For instance, update the relevant DOFs of a discontinuous P1 space
  // (whose DOFs are element-based) after a pressure point
  // has marked a vertex essential.
  //
  // For now we only synchronize vertices and elements DOF (e.g. to apply pressure point (0D) BC),
  // additional checks can be added when necessary 
  // (e.g., edges and elements DOF for Dirichlet BC on 1D boundaries).
  for(feSpace *fS : spaces) {
    fS->synchronizeCodeOfEssentialDOF();
    // fS->synchronizeNumberingOfEssentialDOF();
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
  _nInc = globalNum;
  for(int i = 0; i < _nFields; ++i) {
    globalNum = _numberings[_fieldIDs[i]]->numberEssential(globalNum);
  }
  _nDofs = globalNum;

  // Now synchronize the actual DOF numbers of discontinuous spaces
  // Keep track of the number of DOF whose number has been overriden,
  // as the number of all other DOF with higher number must be decremented
  // to avoid holes in the global numbering
  int numModifiedDOF = 0;
  for(feSpace *fS : spaces) {
    fS->synchronizeNumberingOfEssentialDOF(numModifiedDOF);
  }
  if(numModifiedDOF > 0) {
    for(int i = 0; i < _nFields; ++i) {
      int highestUnmodifiedDOF = _nDofs - numModifiedDOF - 1;
      for(auto &dof : _numberings[_fieldIDs[i]]->_numberingVertices) {
        if(dof > highestUnmodifiedDOF) dof -= numModifiedDOF;
      }
      for(auto &dof : _numberings[_fieldIDs[i]]->_numberingElements) {
        if(dof > highestUnmodifiedDOF) dof -= numModifiedDOF;
      }
      for(auto &dof : _numberings[_fieldIDs[i]]->_numberingEdges) {
        if(dof > highestUnmodifiedDOF) dof -= numModifiedDOF;
      }
      for(auto &dof : _numberings[_fieldIDs[i]]->_numberingFaces) {
        if(dof > highestUnmodifiedDOF) dof -= numModifiedDOF;
      }
    }
    _nDofs -= numModifiedDOF;
  }

  // Handle periodicity
  for(const auto &fieldID : _fieldIDs)
  {
    _numberings[fieldID]->applyPeriodicity(mesh, spaces);

    for(const auto &pair : _numberings[fieldID]->PeriodicDOF())
    {
      _periodicDOF.insert(pair);
    }
  }

  // Summarize all the degrees of freedom assigned for each field in the _allFieldDOF vector
  for(int i = 0; i < _nFields; ++i) {
    _numberings[_fieldIDs[i]]->compactFieldDOF();
  }

  for(feSpace *fS : spaces) {
    fS->setAsNumbered();
  }

  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "DEGREES OF FREEDOM:");
  for(int i = 0; i < _nFields; ++i) {
    feInfoCond(FE_VERBOSE > 0, "\t\tInfo for field \"%s\":", _fieldIDs[i].data());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of DOF: %d", _numberings[_fieldIDs[i]]->getNbDOFs());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of unknown DOF: %d",
               _numberings[_fieldIDs[i]]->getUnknownDOF().size());
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of essential DOF: %d",
               _numberings[_fieldIDs[i]]->getEssentialDOF().size());
  }
}

feMetaNumber::~feMetaNumber()
{
  // Delete all numberings
  for(auto it = _numberings.begin(); it != _numberings.end();) {
    delete it->second;
    it = _numberings.erase(it);
  }
}

feStatus feMetaNumber::exportNumberingVertices(feMesh *mesh,
  std::string fileName,
  bool posFormat)
{
  FILE *file = fopen(fileName.data(), "w");
  if(file == nullptr) {
    return feErrorMsg(FE_STATUS_WRITE_ERROR, "Could not export field numbering to file %s",
                      fileName.data());
  }

  if(posFormat)
  {
    fprintf(file, "View \"%s\"{\n", fileName.data());
  }

  for(auto pair : _numberings) {
    pair.second->exportNumberingVertices(mesh, file, posFormat);
  }

  if(posFormat)
  {
    fprintf(file, "};");
  }

  fclose(file);
  return FE_STATUS_OK;
}

void feMetaNumber::printFields(std::ostream &os)
{
  for(auto s : _fieldIDs) os << s << std::endl;
}

void feMetaNumber::printNumberings(std::ostream &os)
{
  for(int i = 0; i < _nFields; ++i) {
    os << "Field " << _fieldIDs[i] << " - Numbering of vertex degrees of freedom:" << std::endl;
    _numberings[_fieldIDs[i]]->printNumberingVertices(os);
    os << "Field " << _fieldIDs[i] << " - Numbering of element degrees of freedom:" << std::endl;
    _numberings[_fieldIDs[i]]->printNumberingElements(os);
    os << "Field " << _fieldIDs[i] << " - Numbering of edge degrees of freedom:" << std::endl;
    _numberings[_fieldIDs[i]]->printNumberingEdges(os);
  }
}

void feMetaNumber::printCodes(std::ostream &os)
{ 
  // Attention : ordre dépend du nom des variables (ordre du map.first)
  for(auto const &secondIsNumber : _numberings) {
    os << "Field " << secondIsNumber.first << " - vertices :" << std::endl;
    secondIsNumber.second->printCodeVertices(os);
    os << "Field " << secondIsNumber.first << " - elements :" << std::endl;
    secondIsNumber.second->printCodeElements(os);
    os << "Field " << secondIsNumber.first << " - edges :" << std::endl;
    secondIsNumber.second->printCodeEdges(os);
  }
}