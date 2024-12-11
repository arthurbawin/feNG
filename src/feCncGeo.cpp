#include "feNG.h"
#include "feCncGeo.h"
#include "feSpace.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "SPoint2.h"

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
    case geometricInterpolant::TETP1:
      return "TETP1";
    case geometricInterpolant::TETP2:
      return "TETP2";
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
    _connecVertices(connecVertices),  _connecEdges(connecEdges), _connecFaces(connecFaces),
    _connecElem(connecElem), _geometricInterpolant(space)
{
  if(connecElem.size() == 0) _connecElem.resize(nElements);
  if(connecEdges.size() == 0) _connecEdges.resize(nElements * nEdgesPerElement);

  _elementsVolume.resize(_nElements, 0.);
  _minimumScaledJacobianControlCoeffs.resize(nElements, 0.);

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
  colorElements(3);
  printColoringStatistics();

  // _mycoloring = new feColoring(1, _nVerticesPerElm, _connecVertices, _connecElem);
}

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
// #if defined(FENG_DEBUG)
  if((size_t) (_nVerticesPerElm * numElem + iVertex) >= _connecVertices.size())
    feErrorMsg(FE_STATUS_ERROR,
               "Out of bounds: accessing entry %d in _connecVertices"
               " of size %u",
               _nVerticesPerElm * numElem + iVertex, _connecVertices.size());
// #endif
  return _connecVertices[_nVerticesPerElm * numElem + iVertex];
}

void feCncGeo::setVertexConnectivity(const int numElem, const int iVertex, const int val)
{
#if defined(FENG_DEBUG)
  if((size_t) _(nVerticesPerElm * numElem + iVertex) >= _connecVertices.size())
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
  if((size_t) numElem >= _connecElem.size())
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
  if((size_t) numElem >= _connecElem.size())
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
  if((size_t) iEdge >= _connecEdgesOnly.size())
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
  if((size_t) (_nEdgesPerElem * numElem + iEdge) >= _connecEdges.size())
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

feStatus feCncGeo::computeJacobians(const bool ignoreNegativeJacobianWarning)
{
  int nQuad = _geometricInterpolant->getNumQuadPoints();
  std::vector<double> &wQuad = _geometricInterpolant->getQuadratureWeights();
  _J.resize(_nElements * nQuad);

  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  // FILE *myfile = fopen("jacobians.pos", "w");
  // fprintf(myfile, "View \"jacobians\"{\n");

  bool atLeastOneNegative = false;

  switch(_dim) {
    case 0:
      for(int iElm = 0; iElm < _nElements; ++iElm) {
        for(int k = 0; k < nQuad; ++k) {
          _J[nQuad * iElm + k] = 1.;
        }
        _elementsVolume[iElm] = 1.;
      }
      break;

    case 1: {
      std::vector<double> dxdr(3, 0.);

      for(int iElm = 0; iElm < _nElements; ++iElm) {
        _mesh->getCoord(_tag, iElm, geoCoord);
        _elementsVolume[iElm] = 0.;
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _J[nQuad * iElm + k] = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
          _elementsVolume[iElm] += wQuad[k] * _J[nQuad * iElm + k];
        }
      }
      break;
    }

    case 2: {
      std::vector<double> dxdr(3, 0.); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.); // [dx/ds, dy/ds, dz/ds]
      for(int iElm = 0; iElm < _nElements; ++iElm)
      {
        _mesh->getCoord(_tag, iElm, geoCoord);
        _elementsVolume[iElm] = 0.;

        double jMin = 1e22;
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
          _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
          _elementsVolume[iElm] += wQuad[k] * _J[nQuad * iElm + k];
          jMin = fmin(jMin, _J[nQuad * iElm + k]);

          if(_J[nQuad * iElm + k] <= 0) {
            atLeastOneNegative = true;
            if(_nVerticesPerElm == 2) {
              feWarning("Negative or zero jacobian = %+-12.12e on elm %d with coordinates "
                         "(%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e)\n",
                         _J[nQuad * iElm + k], iElm,
                         geoCoord[3 * 0 + 0], geoCoord[3 * 0 + 1],
                         geoCoord[3 * 1 + 0], geoCoord[3 * 1 + 1]);
            } else if(_nVerticesPerElm == 3) {
              feWarning("Negative or zero jacobian = %+-12.12e on elm %d with coordinates "
                       "(%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e)\n",
                       _J[nQuad * iElm + k], iElm,
                       geoCoord[3 * 0 + 0], geoCoord[3 * 0 + 1],
                       geoCoord[3 * 1 + 0], geoCoord[3 * 1 + 1],
                       geoCoord[3 * 2 + 0], geoCoord[3 * 2 + 1]);
            } else if(_nVerticesPerElm == 6) {
              feWarning("Negative or zero jacobian = %+-12.12e on elm %d with coordinates "
                       "(%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e) - "
                       "(%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e) - (%+-1.4e - %+-1.4e)\n",
                       _J[nQuad * iElm + k], iElm, 
                       geoCoord[3 * 0 + 0], geoCoord[3 * 0 + 1],
                       geoCoord[3 * 1 + 0], geoCoord[3 * 1 + 1],
                       geoCoord[3 * 2 + 0], geoCoord[3 * 2 + 1],
                       geoCoord[3 * 3 + 0], geoCoord[3 * 3 + 1],
                       geoCoord[3 * 4 + 0], geoCoord[3 * 4 + 1],
                       geoCoord[3 * 5 + 0], geoCoord[3 * 5 + 1]);
            }
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

  // fprintf(myfile, "};\n"); fclose(myfile);

  if(atLeastOneNegative && !ignoreNegativeJacobianWarning) {
    return feErrorMsg(FE_STATUS_ERROR, "Negative or zero jacobian on at least one element )-:");
  }

  return FE_STATUS_OK;
}

feStatus feCncGeo::recomputeElementJacobian(const int iElm)
{
  int nQuad = _geometricInterpolant->getNumQuadPoints();
  std::vector<double> &wQuad = _geometricInterpolant->getQuadratureWeights();
  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  bool atLeastOneNegative = false;

  switch(_dim) {
    case 0:
      for(int i = 0; i < _nElements; ++i) {
        for(int k = 0; k < nQuad; ++k) {
          _J[nQuad * i + k] = 1.0;
        }
      }
      break;

    case 1: {
      std::vector<double> dxdr(3, 0.0);

      for(int i = 0; i < _nElements; ++i) {
        _mesh->getCoord(_tag, i, geoCoord);
        _elementsVolume[i] = 0.;
        for(int k = 0; k < nQuad; ++k) {
          _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
          _J[nQuad * i + k] = sqrt(dxdr[0] * dxdr[0] + dxdr[1] * dxdr[1]);
          _elementsVolume[i] += wQuad[k] * _J[nQuad * i + k];
        }
      }
      break;
    }

    case 2: {
      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]

      _mesh->getCoord(_tag, iElm, geoCoord);
      _elementsVolume[iElm] = 0.;

      double jMin = 1e22;
      for(int k = 0; k < nQuad; ++k) {
        _geometricInterpolant->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
        _geometricInterpolant->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
        _J[nQuad * iElm + k] = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
        _elementsVolume[iElm] += wQuad[k] * _J[nQuad * iElm + k];
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

feStatus feCncGeo::computeNormalVectors(std::vector<double> &normalVectors) const
{
  int nQuad = _geometricInterpolant->getNumQuadPoints();
  std::vector<double> &xQuad = _geometricInterpolant->getRQuadraturePoints();
  normalVectors.resize(3 * _nElements * nQuad);

  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  if(_dim != 1) {
    return feErrorMsg(FE_STATUS_ERROR, "Only computing normal vectors for 1D connectivities for now.");
  }

  for(int iElm = 0; iElm < _nElements; ++iElm) {
    _mesh->getCoord(_tag, iElm, geoCoord);
    for(int k = 0; k < nQuad; ++k) {
      double tangent[2] = {1., 0.};

      if(_interpolant == geometricInterpolant::LINEP1) {
        tangent[0] = geoCoord[3] - geoCoord[0];
        tangent[1] = geoCoord[4] - geoCoord[1];
      }
      if(_interpolant == geometricInterpolant::LINEP2) {
        // Map quad r-coordinate from [-1,1] to [0,1] 
        double t = (xQuad[k] + 1.) / 2.;
        tangent[0] = geoCoord[0] * (4.*t - 3.) + geoCoord[6] * (4. - 8.*t) + geoCoord[3] * (4.*t - 1.);
        tangent[1] = geoCoord[1] * (4.*t - 3.) + geoCoord[7] * (4. - 8.*t) + geoCoord[4] * (4.*t - 1.);
      }

      double norm = sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
      normalVectors[3 * nQuad * iElm + 3 * k + 0] = -tangent[1] / norm;
      normalVectors[3 * nQuad * iElm + 3 * k + 1] =  tangent[0] / norm;
      normalVectors[3 * nQuad * iElm + 3 * k + 2] = 0.;
    }
  }

  return FE_STATUS_OK;
}

// Lagrange points are given as {L200, L020, L002, L110, L011, L101}
void getBezierControlPoints(const SPoint2 L[6], SPoint2 P[6])
{
  P[0] = L[0]; // P_200
  P[1] = L[1]; // P_020
  P[2] = L[2]; // P_002
  P[3] = (L[0] * (-1.) + L[3] * 4. - L[1]) * 0.5; // P_110
  P[4] = (L[1] * (-1.) + L[4] * 4. - L[2]) * 0.5; // P_011
  P[5] = (L[2] * (-1.) + L[5] * 4. - L[0]) * 0.5; // P_101
}

// Compute determinant of the 2x2 matrix |P0 P1|,
// points are column vectors.
inline double computeDeterminant(const SPoint2 &P0, const SPoint2 &P1)
{
  return (P0[0] * P1[1]) - (P1[0] * P0[1]);
}

double controlCoefficientN200(const SPoint2 &P200, const SPoint2 &P011,
                              const SPoint2 &P020, const SPoint2 &P101,
                              const SPoint2 &P002, const SPoint2 &P110)
{
  UNUSED(P011, P020, P002);
  return 4. * computeDeterminant(P200 - P110, P200 - P101);
}

double controlCoefficientN110(const SPoint2 &P200, const SPoint2 &P011,
                              const SPoint2 &P020, const SPoint2 &P101,
                              const SPoint2 &P002, const SPoint2 &P110)
{
  UNUSED(P002);
  return 2. * computeDeterminant(P200 - P101, P020 - P011) + 2. * computeDeterminant(P110 - P011, P110 - P101);
}

// Bezier points are given as {P200, P020, P002, P110, P011, P101}
void getJacobianControlCoefficients(const SPoint2 P[6], double N[6])
{
  const SPoint2 &P200 = P[0];
  const SPoint2 &P020 = P[1];
  const SPoint2 &P002 = P[2];
  const SPoint2 &P110 = P[3];
  const SPoint2 &P011 = P[4];
  const SPoint2 &P101 = P[5];
  N[0] = controlCoefficientN200(P200, P011, P020, P101, P002, P110);
  N[1] = controlCoefficientN200(P020, P101, P002, P110, P200, P011);
  N[2] = controlCoefficientN200(P002, P110, P200, P011, P020, P101);
  N[3] = controlCoefficientN110(P200, P011, P020, P101, P002, P110);
  N[4] = controlCoefficientN110(P020, P101, P002, P110, P200, P011);
  N[5] = controlCoefficientN110(P002, P110, P200, P011, P020, P101);
}

// uvw are the barycentric coordinates: in 2D triangle:
// uvw[0] = 1. - xsi - eta;
// uvw[1] = xsi;
// uvw[2] = eta;
void getP2BernsteinBasis(const double uvw[3], double basis[6])
{
  double u = uvw[0];
  double v = uvw[1];
  double w = uvw[2];
  basis[0] = u*u;
  basis[1] = 2.*u*v;
  basis[2] = v*v;
  basis[3] = 2.*v*w;
  basis[4] = w*w;
  basis[5] = 2.*w*u;
}

feStatus feCncGeo::computeMinimumScaledJacobianControlCoefficients()
{
  if(_dim != 2 || _nVerticesPerElm != 6) {
    return feErrorMsg(FE_STATUS_ERROR, "Only computing scaled jacobian control coefficients for P2 triangles.");
  }

  std::vector<double> geoCoord(3 * 6);
  SPoint2 lagrangePoints[6], bezierPoints[6];
  double N[6];

  for(int iElm = 0; iElm < _nElements; ++iElm)
  {
    _mesh->getCoord(_tag, iElm, geoCoord);

    // Get Lagrange and Bezier control points
    for(int i = 0; i < 6; ++i) {
      lagrangePoints[i] = SPoint2(geoCoord[i * 3 + 0], geoCoord[i * 3 + 1]);
    }
    getBezierControlPoints(lagrangePoints, bezierPoints);

    // Get control coefficients
    getJacobianControlCoefficients(bezierPoints, N);

    double minControlCoeff = 1e22;
    for(int i = 0; i < 6; ++i) {
      minControlCoeff = fmin(minControlCoeff, N[i]);
    }

    _minimumScaledJacobianControlCoeffs[iElm] = minControlCoeff/(2. * _elementsVolume[iElm]);
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
  tic();

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

  feInfoCond(FE_VERBOSE > 0, "\t\tDone in %f s", toc());
}

void feCncGeo::printColoringStatistics()
{
  feInfoCond(FE_VERBOSE > 0, "\t\tColoring statistics:");
  for(size_t i = 0; i < _nbElmPerColor.size(); ++i) {
    feInfoCond(FE_VERBOSE > 0, "\t\t\tNumber of elements in color %d: %d", i, _nbElmPerColor[i]);
  }
}

// From built-in routines to remove
void feCncGeo::printColoring(std::string fileName)
{
  FILE *myfile = fopen(fileName.data(), "w");
  fprintf(myfile, "View \"%s\"{\n", fileName.data());
  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  for(size_t iElm = 0; iElm < _elmToColor.size(); ++iElm) {
    _mesh->getCoord(_tag, iElm, geoCoord);
    writeElementToPOS(myfile, geoCoord, (double) _elmToColor[iElm]);
  }
  fprintf(myfile, "};\n"); fclose(myfile);
}

// From class feColoring
void feCncGeo::printColoring2(std::string fileName)
{
  FILE *myfile = fopen(fileName.data(), "w");
  fprintf(myfile, "View \"%s\"{\n", fileName.data());
  std::vector<double> geoCoord(3 * _nVerticesPerElm);

  // const std::vector<int> &_elemToColor = _mycoloring->getColorElm();

  for(size_t iElm = 0; iElm < _elmToColor.size(); ++iElm) {
    _mesh->getCoord(_tag, iElm, geoCoord);
    writeElementToPOS(myfile, geoCoord, (double) _elmToColor[iElm]);
  }
  fprintf(myfile, "};\n"); fclose(myfile);
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