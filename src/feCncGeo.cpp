#include "feNG.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feMesh.h"
#include "feQuadrature.h"

extern int FE_VERBOSE;

std::string toString(geometryType t)
{
  switch(t) {
    case geometryType::POINT:
      return "POINT";
    case geometryType::LINE:
      return "LINE";
    case geometryType::TRI:
      return "TRI";
    case geometryType::QUAD:
      return "QUAD";
    case geometryType::TET:
      return "TET";
    case geometryType::HEX:
      return "HEX";
    default:
      return "Unknown geometry";
  }
}

std::string toString(geometricInterpolant t)
{
  switch(t) {
    case geometricInterpolant::NONE:
      return "NONE";
    case geometricInterpolant::POINTP0:
      return "POINTP0";
    case geometricInterpolant::LINEP1:
      return "LINEP1";
    case geometricInterpolant::LINEP2:
      return "LINEP2";
    case geometricInterpolant::LINEP3:
      return "LINEP3";
    case geometricInterpolant::TRIP1:
      return "TRIP1";
    case geometricInterpolant::TRIP2:
      return "TRIP2";
    case geometricInterpolant::TRIP3:
      return "TRIP3";
    default:
      return "Unknown geometric interpolant";
  }
}

int getGeometricInterpolantDegree(geometricInterpolant t)
{
  switch(t) {
    case geometricInterpolant::NONE:
      return -1;
    case geometricInterpolant::POINTP0:
      return 0;
    case geometricInterpolant::LINEP1:
      return 1;
    case geometricInterpolant::LINEP2:
      return 2;
    case geometricInterpolant::LINEP3:
      return 3;
    case geometricInterpolant::TRIP1:
      return 1;
    case geometricInterpolant::TRIP2:
      return 2;
    case geometricInterpolant::TRIP3:
      return 3;
    default:
      return -1;
  }
}

feCncGeo::feCncGeo(const int tag, const int dimension, const int nVerticesPerElement,
                   const int nElements, const int nEdgesPerElement, const std::string &ID,
                   const geometryType geometry, const geometricInterpolant interpolant,
                   feSpace *space, std::vector<int> connecVertices, std::vector<int> connecElem,
                   std::vector<int> connecEdges, std::vector<int> connecFaces)
  : _ID(ID), _tag(tag), _dim(dimension), _geometry(geometry), _interpolant(interpolant),
    _nVerticesPerElm(nVerticesPerElement), _nElements(nElements), _nEdgesPerElem(nEdgesPerElement),
    _connecVertices(connecVertices), _connecElem(connecElem), _connecEdges(connecEdges),
    _connecFaces(connecFaces), _geometricInterpolant(space)
{
  if(connecElem.size() == 0) _connecElem.resize(nElements);
  if(connecEdges.size() == 0) _connecEdges.resize(nElements * nEdgesPerElement);

  std::sort(connecVertices.begin(), connecVertices.end());
  _nVertices = std::unique(connecVertices.begin(), connecVertices.end()) - connecVertices.begin();

  // Create the connectivity of unique vertices
  // (_connecVertices has size nElm * nVerticesperElm)
  _connecVerticesOnly = connecVertices;
  std::sort(_connecVerticesOnly.begin(), _connecVerticesOnly.end());
  auto last = std::unique(_connecVerticesOnly.begin(), _connecVerticesOnly.end());
  _connecVerticesOnly.erase(last, _connecVerticesOnly.end());

  // Create the connectivity of unique edges
  _connecEdgesOnly = connecEdges;
  for(size_t i = 0; i < _connecEdgesOnly.size(); ++i)
    _connecEdgesOnly[i] = fabs(_connecEdgesOnly[i]);
  std::sort(_connecEdgesOnly.begin(), _connecEdgesOnly.end());
  last = std::unique(_connecEdgesOnly.begin(), _connecEdgesOnly.end());
  _connecEdgesOnly.erase(last, _connecEdgesOnly.end());
  _nEdges = _connecEdgesOnly.size();

  // Color the elements for partitioning
  colorElements(1);
};

int feCncGeo::getUniqueVertexConnectivity(const int iVertex) const
{
#if defined(FENG_DEBUG)
  if(iVertex >= _connecVerticesOnly.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecVerticesOnly"
               " of size %u",
               iVertex, _connecVerticesOnly.size());
#endif
  return _connecVerticesOnly[iVertex];
}

int feCncGeo::getVertexConnectivity(const int iVertex) const
{
#if defined(FENG_DEBUG)
  if(iVertex >= _connecVertices.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecVertices"
               " of size %u",
               iVertex, _connecVertices.size());
#endif
  return _connecVertices[iVertex];
}

int feCncGeo::getVertexConnectivity(const int numElem, const int iVertex) const
{
#if defined(FENG_DEBUG)
  if(_nVerticesPerElm * numElem + iVertex >= _connecVertices.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecVertices"
               " of size %u",
               _nVerticesPerElm * numElem + iVertex, _connecVertices.size());
#endif
  return _connecVertices[_nVerticesPerElm * numElem + iVertex];
}

void feCncGeo::setVertexConnectivity(const int numElem, const int iVertex, const int val)
{
#if defined(FENG_DEBUG)
  if(_nVerticesPerElm * numElem + iVertex >= _connecVertices.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecVertices"
               " of size %u",
               _nVerticesPerElm * numElem + iVertex, _connecVertices.size());
#endif
  _connecVertices[_nVerticesPerElm * numElem + iVertex] = val;
}

int feCncGeo::getElementConnectivity(const int numElem) const
{
#if defined(FENG_DEBUG)
  if(numElem >= _connecElem.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecElem"
               " of size %u",
               numElem, _connecElem.size());
#endif
  return _connecElem[numElem];
}

void feCncGeo::setElementConnectivity(const int numElem, const int val)
{
#if defined(FENG_DEBUG)
  if(numElem >= _connecElem.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecElem"
               " of size %u",
               numElem, _connecElem.size());
#endif
  _connecElem[numElem] = val;
}

int feCncGeo::getUniqueEdgeConnectivity(const int iEdge) const
{
#if defined(FENG_DEBUG)
  if(iEdge >= _connecEdgesOnly.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecEdgesOnly"
               " of size %u",
               iEdge, _connecEdgesOnly.size());
#endif
  return _connecEdgesOnly[iEdge];
}

int feCncGeo::getEdgeConnectivity(const int numElem, const int iEdge) const
{
#if defined(FENG_DEBUG)
  if(_nEdgesPerElem * numElem + iEdge >= _connecEdges.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecEdges"
               " of size %u",
               _nEdgesPerElem * numElem + iEdge, _connecEdges.size());
#endif
  return _connecEdges[_nEdgesPerElem * numElem + iEdge];
}

void feCncGeo::setEdgeConnectivity(const int numElem, const int iEdge, const int val)
{
#if defined(FENG_DEBUG)
  if(_nEdgesPerElem * numElem + iEdge >= _connecEdges.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecEdges"
               " of size %u",
               _nEdgesPerElem * numElem + iEdge, _connecEdges.size());
#endif
  _connecEdges[_nEdgesPerElem * numElem + iEdge] = val;
}

feStatus feCncGeo::setQuadratureRule(feQuadrature *rule)
{
  return _geometricInterpolant->setQuadratureRule(rule);
}

feStatus feCncGeo::computeJacobians()
{
  int nQuad = _geometricInterpolant->getNumQuadPoints();
  _J.resize(_nElements * nQuad);

  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  // FILE *myfile = fopen("jacobians.pos", "w");
  // fprintf(myfile, "View \"jacobians\"{\n");

  bool atLeastOneNegative = false;

  switch(_dim) {
    case 0:
      for(int iElm = 0; iElm < _nElements; ++iElm) {
        for(int k = 0; k < nQuad; ++k) {
          _J[nQuad * iElm + k] = 1.0;
        }
      }
      break;

    case 1: {
      std::vector<double> dxdr(3, 0.0);

      for(int iElm = 0; iElm < _nElements; ++iElm) {
        _mesh->getCoord(_tag, iElm, geoCoord);
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _J[nQuad * iElm + k] = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
        }
      }
      break;
    }

    case 2: {
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      for(int iElm = 0; iElm < _nElements; ++iElm)
      {
        _mesh->getCoord(_tag, iElm, geoCoord);

        double jMin = 1e22;
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
          _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
          jMin = fmin(jMin, _J[nQuad * iElm + k]);

          if(_J[nQuad * iElm + k] <= 0) {
            atLeastOneNegative = true;
            feWarning("Negative or zero jacobian = %+-12.12e on elm %d with coordinates "
                       "(%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e)\n",
                       _J[nQuad * iElm + k], iElm, geoCoord[3 * 0 + 0], geoCoord[3 * 0 + 1],
                       geoCoord[3 * 1 + 0], geoCoord[3 * 1 + 1], geoCoord[3 * 2 + 0],
                       geoCoord[3 * 2 + 1]);
          }
        }

        // Plot jacobian
        // writeElementToPOS(myfile, geoCoord, jMin);
      }
      break;
    }

    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Element jacobian not implemented "
                        "for elements with dim = %d.\n",
                        _dim);
  }

  if(atLeastOneNegative){
    return feErrorMsg(FE_STATUS_ERROR, "Negative or zero jacobian on at least one element )-:");
  }

  // fprintf(myfile, "};\n"); fclose(myfile);

  return FE_STATUS_OK;
}

feStatus feCncGeo::recomputeElementJacobian(const int iElm)
{
  int nQuad = _geometricInterpolant->getNumQuadPoints();
  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  bool atLeastOneNegative = false;

  switch(_dim) {
    case 0:
      for(int iElm = 0; iElm < _nElements; ++iElm) {
        for(int k = 0; k < nQuad; ++k) {
          _J[nQuad * iElm + k] = 1.0;
        }
      }
      break;

    case 1: {
      std::vector<double> dxdr(3, 0.0);

      for(int iElm = 0; iElm < _nElements; ++iElm) {
        _mesh->getCoord(_tag, iElm, geoCoord);
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _J[nQuad * iElm + k] = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
        }
      }
      break;
    }

    case 2: {
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

      _mesh->getCoord(_tag, iElm, geoCoord);

      double jMin = 1e22;
      for(int k = 0; k < nQuad; ++k) {
        _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
        jMin = fmin(jMin, _J[nQuad * iElm + k]);

        if(_J[nQuad * iElm + k] <= 0) {
          atLeastOneNegative = true;
        }
      }

      break;
    }

    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Element jacobian not implemented "
                        "for elements with dim = %d.\n",
                        _dim);
  }

  if(atLeastOneNegative){
    return FE_STATUS_FAILED;
  }

  return FE_STATUS_OK;
}

void feCncGeo::computeElementTransformation(std::vector<double> &elementCoord, const int iQuadNode,
                                            const double jac, ElementTransformation &T) const
{
  T.jac = jac;
  if(_dim == 0) {
    T.drdx[0] = 1.;
  } else if(_dim == 1) {
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(elementCoord, iQuadNode,
                                                                        T.dxdr);
    T.drdx[0] = 1. / T.dxdr[0];
  } else if(_dim == 2) {
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(elementCoord, iQuadNode,
                                                                        T.dxdr);
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(elementCoord, iQuadNode,
                                                                        T.dxds);
    T.drdx[0] = T.dxds[1] / jac;
    T.drdx[1] = -T.dxdr[1] / jac;
    T.drdy[0] = -T.dxds[0] / jac;
    T.drdy[1] = T.dxdr[0] / jac;
  } else if(_dim == 3) {
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(elementCoord, iQuadNode,
                                                                        T.dxdr);
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(elementCoord, iQuadNode,
                                                                        T.dxds);
    _geometricInterpolant->interpolateVectorFieldAtQuadNode_tDerivative(elementCoord, iQuadNode,
                                                                        T.dxdt);
    /* Write inverse of 3x3 matrix (-: */
    feErrorMsg(FE_STATUS_ERROR, "Inverse of element transformation not implemented in 3D");
    exit(-1);
  }
}

static int colorChoice(std::vector<bool> availableColor, std::vector<int> &nbElmPerColor,
                       int &nbColor, std::vector<std::vector<int> > &listElmPerColor)
{
  int minColor = *max_element(nbElmPerColor.begin(), nbElmPerColor.end());
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true) {
      if(minColor > nbElmPerColor[i]) {
        minColor = nbElmPerColor[i];
      }
    }
  }

  int activColor = -1;
  for(int i = 0; i < nbColor; i++) {
    if(availableColor[i] == true && nbElmPerColor[i] == minColor) {
      activColor = i;
    }
  }

  if(activColor == -1) {
    activColor = nbColor;
    nbColor = nbColor + 1;
    nbElmPerColor.push_back(0);
    listElmPerColor.push_back(std::vector<int>(0, 0));
  }
  return activColor;
}

void feCncGeo::colorElements(int coloringAlgorithm)
{
  feInfoCond(FE_VERBOSE > 0, "\t\tColoring connectivity %s...", _ID.c_str());
  // Create the node patch :
  int size = *std::max_element(_connecVertices.begin(), _connecVertices.end()) + 1;
  _listElmPerNode.resize(size);
  _nbElmPerNode.resize(size);
  for(int i = 0; i < _nElements; ++i) {
    for(int j = 0; j < _nVerticesPerElm; ++j) {
      int nds = _connecVertices[i * _nVerticesPerElm + j];
      _nbElmPerNode[nds] += 1;
      _listElmPerNode[nds].push_back(i);
    }
  }

  switch(coloringAlgorithm) {
    case 1: // Returning a non homogeneous distribution of elements per color
    {
      _nbColor = 0;
      _elmToColor.resize(_nElements, -1);

      bool noColor = false;
      while(!noColor) {
        for(int k = 0; k < _nElements; ++k) {
          if(_elmToColor[k] == -1) {
            for(int i = 0; i < _nVerticesPerElm; ++i) {
              int s = _connecVertices[k * _nVerticesPerElm + i];
              for(int j = 0; j < _nbElmPerNode[s]; ++j) {
                int numElm = _listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elmToColor[numElm] < 0) {
                    _elmToColor[numElm] = -2;
                  }
                }
              }
            }
          }
        }

        _nbElmPerColor.push_back(0);
        _listElmPerColor.push_back(std::vector<int>(0, 0));

        noColor = true;
        for(int k = 0; k < _nElements; ++k) {
          if(_elmToColor[k] == -1) {
            _elmToColor[k] = _nbColor;
            _nbElmPerColor[_nbColor] = _nbElmPerColor[_nbColor] + 1;
            _listElmPerColor[_nbColor].push_back(k);
          } else if(_elmToColor[k] == -2) {
            _elmToColor[k] = -1;
            noColor = false;
          }
        }
        _nbColor = _nbColor + 1;

      } // while
      break;
    }

    case 2: // Returning a homogeneous distribution of elements per color
    {
      // Create the Elements patch :
      int sizeElm = *std::max_element(_connecElem.begin(), _connecElem.end()) + 1;
      _nbElmPerElm.resize(sizeElm);
      _listElmPerElm.resize(sizeElm);
      std::vector<std::set<int> > listElmPerElm(sizeElm);
      for(int i = 0; i < _nElements; ++i) {
        int elm = _connecElem[i];
        for(int j = 0; j < _nVerticesPerElm; ++j) {
          int nds = _connecVertices[i * _nVerticesPerElm + j];
          for(int ngb : _listElmPerNode[nds]) {
            if(ngb != i) {
              listElmPerElm[elm].insert(ngb);
            }
          }
        }
      }
      // Convert sets to vectors
      for(int i = 0; i < _nElements; ++i) {
        int elm = _connecElem[i];
        std::vector<int> v(listElmPerElm[elm].begin(), listElmPerElm[elm].end());
        _nbElmPerElm[elm] = v.size();
        _listElmPerElm[elm] = v;
      }

      // Finding the Elm with the most neighbours to asign nb Colors
      _nbColor = *std::max_element(_nbElmPerElm.begin(), _nbElmPerElm.end()) + 1;

      _elmToColor.resize(_nElements, -1);
      _nbElmPerColor.resize(_nbColor);
      _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));
      for(int iElm = 0; iElm < _nElements; iElm++) {
        if(_elmToColor[iElm] == -1) {
          int elm = _connecElem[iElm];
          // feInfo("%d",elm);
          std::vector<int> iPatchElm = _listElmPerElm[elm];
          std::vector<bool> availableColor(_nbColor, true);
          for(int i = 0; i < _nbElmPerElm[elm]; i++) {
            if(_elmToColor[iPatchElm[i]] != -1 &&
               availableColor[_elmToColor[iPatchElm[i]]] == true) {
              availableColor[_elmToColor[iPatchElm[i]]] = false;
            }
          }

          int activColor = colorChoice(availableColor, _nbElmPerColor, _nbColor, _listElmPerColor);
          _elmToColor[iElm] = activColor;
          _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
          _listElmPerColor[activColor].push_back(iElm);
        }
      }
      break;
    }

    case 3: // Idem 2 but not using Patch Elm
    {
      // Finding the Node with the most neighbours to asign initial nb Colors
      _nbColor = *std::max_element(_nbElmPerNode.begin(), _nbElmPerNode.end());
      // for(auto pair : _nbElmPerNode) {
      //   if(pair.second > _nbColor) {
      //     _nbColor = pair.second;
      //   }
      // }

      _elmToColor.resize(_nElements, -1);
      _nbElmPerColor.resize(_nbColor);
      _listElmPerColor.resize(_nbColor, std::vector<int>(0, 0));

      bool noColor = false;
      while(noColor == false) {
        for(int k = 0; k < _nElements; ++k) {
          if(_elmToColor[k] == -1) {
            std::vector<bool> availableColor(_nbColor, true);
            for(int i = 0; i < _nVerticesPerElm; ++i) {
              int s = _connecVertices[k * _nVerticesPerElm + i];
              for(int j = 0; j < _nbElmPerNode[s]; ++j) {
                int numElm = _listElmPerNode[s][j];
                if(numElm != k) {
                  if(_elmToColor[numElm] < 0) {
                    _elmToColor[numElm] = -2;
                  } else if(availableColor[_elmToColor[numElm]] == true) {
                    availableColor[_elmToColor[numElm]] = false;
                  }
                }
              }
            }
            int activColor =
              colorChoice(availableColor, _nbElmPerColor, _nbColor, _listElmPerColor);
            _elmToColor[k] = activColor;
            // feInfo("%d couleur pour elm %d",activColor,k);
            _nbElmPerColor[activColor] = _nbElmPerColor[activColor] + 1;
            _listElmPerColor[activColor].push_back(k);
          }
        }

        noColor = true;
        for(int k = 0; k < _nElements; ++k) {
          if(_elmToColor[k] == -2) {
            _elmToColor[k] = -1;
            noColor = false;
          }
        }
      } // while
      break;
    }
  }

  _coloring.numColors = _nbColor;
  _coloring.elem2Color = _elmToColor;
  _coloring.numElemPerColor = _nbElmPerColor;
  _coloring.elementsInColor = _listElmPerColor;

  feInfoCond(FE_VERBOSE > 0, "\t\tDone");
}

void feCncGeo::writeElementToPOS(FILE *posFile, const std::vector<double> &elementCoord,
                                 const double value) const
{
  switch(_interpolant) {
    case geometricInterpolant::NONE:
      // Do nothing
      break;
    case geometricInterpolant::POINTP0:
      fprintf(posFile, "SP(%g,%g,%g){%g};\n", elementCoord[0], elementCoord[1], elementCoord[2],
              value);
      break;
    case geometricInterpolant::LINEP1:
      fprintf(posFile, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n", elementCoord[0], elementCoord[1],
              elementCoord[2], elementCoord[3], elementCoord[4], elementCoord[5], value, value);
      break;
    case geometricInterpolant::LINEP2:
      fprintf(posFile, "SL2(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n", elementCoord[0],
              elementCoord[1], elementCoord[2], elementCoord[3], elementCoord[4], elementCoord[5],
              elementCoord[6], elementCoord[7], elementCoord[8], value, value, value);
      break;
    case geometricInterpolant::TRIP1:
      fprintf(posFile, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n", elementCoord[0],
              elementCoord[1], elementCoord[2], elementCoord[3], elementCoord[4], elementCoord[5],
              elementCoord[6], elementCoord[7], elementCoord[8], value, value, value);
      break;
    case geometricInterpolant::TRIP2:
      fprintf(posFile,
              "ST2(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g,%g,%g,%g};\n",
              elementCoord[0], elementCoord[1], elementCoord[2], elementCoord[3], elementCoord[4],
              elementCoord[5], elementCoord[6], elementCoord[7], elementCoord[8], elementCoord[9],
              elementCoord[10], elementCoord[11], elementCoord[12], elementCoord[13],
              elementCoord[14], elementCoord[15], elementCoord[16], elementCoord[17], value, value,
              value, value, value, value);
      break;
    default:
      feWarning(
        "Cannot write element to POS file because geometry and/or order is not supported :/");
  }
}