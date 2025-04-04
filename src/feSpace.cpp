#include "feSpace.h"
#include "feMesh.h"
#include "feSolution.h"

extern int FE_VERBOSE;

feStatus createFiniteElementSpace(feSpace *&space, feMesh *mesh, const elementType element,
                                  const int degree, const std::string fieldName,
                                  const std::string cncGeoID, const int degreeQuadrature, void *fct,
                                  const bool useGlobalShapeFunctions)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "FINITE ELEMENT SPACE:");
  feInfoCond(FE_VERBOSE > 0, "\t\tCreating FE space for field \"%s\" on entity \"%s\"",
             fieldName.data(), cncGeoID.data());

  if(mesh == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null mesh pointer.");
  if(fct == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null function pointer.");

  if(mesh->getCncGeoTag(cncGeoID) == -1) {
    return feErrorMsg(FE_STATUS_ERROR, "Geometric connectivity \"%s\" does not exist in the mesh.",
                      cncGeoID.c_str());
  }

  const feCncGeo *cnc = mesh->getCncGeoByName(cncGeoID);
  int dimension = cnc->getDim();
  feInfoCond(FE_VERBOSE > 0, "\t\t\tDimension: %d", dimension);

  geometryType geometry = cnc->getGeometry();

  switch(dimension) {
    case 0:
      if(geometry != geometryType::POINT)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Point");

      if(element == elementType::LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Lagrange");
        space = new feSpace0DP0(mesh, fieldName, cncGeoID, (feFunction *)fct);

      } else if(element == elementType::HERMITE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Hermite");
        space = new feSpace0D_Hermite(mesh, fieldName, cncGeoID, (feFunction *)fct);
      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported finite element.");
      }

      break;

    case 1:
      if(geometry != geometryType::LINE)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Line");

      if(element == elementType::LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Lagrange");
        switch(degree) {
          case 0:
            space = new feSpace1DP0(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 1:
            space = new feSpace1DP1(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 2:
            space = new feSpace1DP2(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 3:
            space = new feSpace1DP3(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 4:
            space = new feSpace1DP4(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No LAGRANGE 1D finite element space implemented for deg > 4.");
        }

      } else if(element == elementType::VECTOR_LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: vector-valued Lagrange");

        switch(degree) {
          case 1:
            space = new feSpaceVecP1<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 2:
            space = new feSpaceVecP2<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 3:
            space = new feSpaceVecP3<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 4:
            space = new feSpaceVecP4<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No vector-valued Lagrange finite element space implemented for deg > 4 or = 0.");
        }
      } else if(element == elementType::LAGRANGE_DISCONTINUOUS) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Discontinuous Lagrange");
        space = new feSpace1D_DG_P1(mesh, fieldName, cncGeoID, (feFunction *)fct);

      } else if(element == elementType::LEGENDRE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Legendre");
        space = new feSpace1D_Legendre(degree, mesh, fieldName, cncGeoID, (feFunction *)fct);

      } else if(element == elementType::CROUZEIX_RAVIART) {
        switch(degree) {
          case 1:
            feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Crouzeix-Raviart");
            space = new feSpace1D_CR0(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No CROUZEIX_RAVIART 1D finite element space implemented for deg > 1.");
        }

      } else if(element == elementType::RAVIART_THOMAS) {
        switch(degree) {
          case 1:
            feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Raviart-Thomas");
            space = new feSpaceTriRT1(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No RAVIART_THOMAS finite element space implemented for deg > 1.");
        }

      } else if(element == elementType::HERMITE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Hermite");
        space = new feSpace1D_H3(mesh, fieldName, cncGeoID, (feFunction *)fct);

      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported finite element.");
      }
      break;

    case 2:
      if(geometry != geometryType::TRI) return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Triangle");

      if(element == elementType::LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Lagrange");

        switch(degree) {
          case 0:
            space = new feSpaceTriP0(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 1:
            space = new feSpaceTriP1(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 2:
            space = new feSpaceTriP2(mesh, fieldName, cncGeoID, (feFunction *)fct,
                                     useGlobalShapeFunctions);
            break;
          case 3:
            space = new feSpaceTriP3(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 4:
            space = new feSpaceTriP4(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No LAGRANGE 2D finite element space implemented on triangles for deg > 4.");
        }

      } else if(element == elementType::LAGRANGE_BUBBLE)
      {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Bubble-enriched Lagrange");
        space = new feSpaceTriP2Bubble(mesh, fieldName, cncGeoID, (feFunction *)fct);

        if(degree != 2) {
          return feErrorMsg(FE_STATUS_ERROR, "Bubble-enriched Lagrange FE space only for quadratic (P2+) elements so far.");
        }

      } else if(element == elementType::LAGRANGE_DISCONTINUOUS) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Discontinuous Lagrange");
        space = new feSpaceTriP1_Discontinuous(mesh, fieldName, cncGeoID, (feFunction *)fct);

        if(degree > 1) {
          return feErrorMsg(FE_STATUS_ERROR, "Discontinuous Lagrange FE space only for linear elements so far.");
        }

      } else if(element == elementType::VECTOR_LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: vector-valued Lagrange");

        switch(degree) {
          case 1:
            space = new feSpaceVecTriP1<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 2:
            space = new feSpaceVecTriP2<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 3:
            space = new feSpaceVecTriP3<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          case 4:
            space = new feSpaceVecTriP4<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No vector-valued finite element space implemented for deg > 4 or = 0.");
        }
      } else if(element == elementType::VECTOR_LAGRANGE_BUBBLE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Bubble-enriched vector-valued Lagrange");

        space = new feSpaceVecTriP2Bubble<2>(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);

        if(degree != 2) {
          return feErrorMsg(FE_STATUS_ERROR, "Bubble-enriched vector Lagrange FE space only for quadratic (P2+) elements so far.");
        }
        
      } else if(element == elementType::CROUZEIX_RAVIART) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: non-conforming Crouzeix-Raviart");

        switch(degree) {
          case 1:
            space = new feSpaceTri_CR1(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          case 2:
            space = new feSpaceTri_CR2(mesh, fieldName, cncGeoID, (feFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No Crouzeix-Raviart finite element space implemented for deg > 2 or = 0.");
        }
      } else if(element == elementType::RAVIART_THOMAS) {
        switch(degree) {
          case 1:
            feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Raviart-Thomas");
            space = new feSpaceTriRT1(mesh, fieldName, cncGeoID, (feVectorFunction *)fct);
            break;
          default:
            return feErrorMsg(
              FE_STATUS_ERROR,
              "No RAVIART_THOMAS finite element space implemented for deg > 1.");
        }

      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported finite element.");
      }
      break;

    case 3:
      if(geometry != geometryType::TET)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Tetrahedron");

      if(element == elementType::LAGRANGE) {
        feInfoCond(FE_VERBOSE > 0, "\t\t\tFinite element: Lagrange");
        space = new feSpaceTetPn(degree, mesh, fieldName, cncGeoID, (feFunction *)fct);
      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported finite element.");
      }
      break;

    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Cannot create a finite element space for dimension > 3 or < 0.");
  }

  feInfoCond(FE_VERBOSE > 0, "\t\t\tPolynomial degree: %d", degree);

  // Set the quadrature rule on this space and the corresponding geometric interpolation space
  feQuadrature rule(degreeQuadrature, mesh->getCncGeoByName(cncGeoID)->getGeometry());

  // FIXME: the FE space of the geometric connectivity is the one of the last created space...
  feStatus s;
  s = space->getCncGeo()->getFeSpace()->setQuadratureRule(&rule);
  if(s != FE_STATUS_OK) { return s; }
  s = space->setQuadratureRule(&rule);
  if(s != FE_STATUS_OK) { return s; }

  return FE_STATUS_OK;
}

feSpace::feSpace(const int dimension, feMesh *mesh, const std::string &fieldID,
                 const std::string &cncGeoID, feFunction *scalarField,
                 feVectorFunction *vectorField, bool useGlobalShapeFunctions)
  : _mesh(mesh), _fieldID(fieldID), _fieldTag(-1), _cncGeoID(cncGeoID),
    _cncGeoTag(-1), _dim(dimension), _nQuad(-1), _nFunctions(0), _useGlobalShapeFunctions(useGlobalShapeFunctions),
    _scalarFct(scalarField), _vectorFct(vectorField)    
{
  // Called with mesh = nullptr when creating the geometric interpolation space
  if(mesh != nullptr) {
    _cnc = mesh->getCncGeoByName(cncGeoID);
    _cncGeoTag = mesh->getCncGeoTag(cncGeoID);
  }
}

feSpace::~feSpace()
{
  // if(_scalarFct) delete _scalarFct;
  // if(_vectorFct) delete _vectorFct;
}

feScalarSpace::feScalarSpace(const int dimension, feMesh *mesh, const std::string &fieldID,
                             const std::string &cncGeoID, feFunction *scalarField,
                             bool useGlobalShapeFunctions)
  : feSpace(dimension, mesh, fieldID, cncGeoID, scalarField, nullptr, useGlobalShapeFunctions)
{
}

int feSpace::getNumElements() { return _mesh->getNumElements(_cncGeoTag); }
int feSpace::getNumVerticesPerElem() { return _mesh->getNumVerticesPerElem(_cncGeoTag); }

feStatus feSpace::setQuadratureRule(feQuadrature *rule)
{
  _nQuad = rule->getNumQuadPoints();
  _wQuad = rule->getWeights();
  _xQuad = rule->getXPoints();
  _yQuad = rule->getYPoints();
  _zQuad = rule->getZPoints();

  // int shapeSize = _nFunctions * _nQuad;
  int shapeSize = _nFunctions * _nComponents * _nQuad;
  _L.resize(shapeSize, 0.0);
  _dLdr.resize(shapeSize, 0.0);
  _dLds.resize(shapeSize, 0.0);
  _dLdt.resize(shapeSize, 0.0);
  _d2Ldr2.resize(shapeSize, 0.0);
  _d2Ldrs.resize(shapeSize, 0.0);
  _d2Lds2.resize(shapeSize, 0.0);
  _d2Ldt2.resize(shapeSize, 0.0);

  // _barycentricCoordinates.resize(3 * _nQuad, 0.);

  /* Reference frame discretization : shape functions are computed once on the reference element,
  then evaluated at the quadrature nodes. */
  for(size_t i = 0; i < _xQuad.size(); ++i) {
    double r[3] = {_xQuad[i], _yQuad[i], _zQuad[i]};
    std::vector<double> l = L(r);
    std::vector<double> dldr = dLdr(r);
    std::vector<double> dlds = dLds(r);
    std::vector<double> dldt = dLdt(r);
    std::vector<double> d2ldr2 = d2Ldr2(r);
    std::vector<double> d2ldrs = d2Ldrs(r);
    std::vector<double> d2lds2 = d2Lds2(r);
    std::vector<double> d2ldt2 = d2Ldt2(r);
    int offset = _nFunctions * _nComponents;
    // int offset = _nFunctions;
    for(int j = 0; j < offset; ++j)
    {
           _L[offset * i + j] = l[j];
        _dLdr[offset * i + j] = dldr[j];
        _dLds[offset * i + j] = dlds[j];
        _dLdt[offset * i + j] = dldt[j];
      _d2Ldr2[offset * i + j] = d2ldr2[j];
      _d2Ldrs[offset * i + j] = d2ldrs[j];
      _d2Lds2[offset * i + j] = d2lds2[j];
      _d2Ldt2[offset * i + j] = d2ldt2[j];
    }

    // _barycentricCoordinates[3 * i + 0] = 1. - r[0] - r[1];
    // _barycentricCoordinates[3 * i + 1] = r[0];
    // _barycentricCoordinates[3 * i + 2] = r[1];
  }

  /* Physical frame discretization : shape functions are computed on the physical element. */
  if(_fieldID != "GEO" && _useGlobalShapeFunctions) {
    feInfo("USING GLOBAL FUNCTIONS on space %s - %s", this->_fieldID.c_str(),
           this->_cncGeoID.c_str());
    
    feCncGeo *cnc = this->getCncGeo();
    int nElm = cnc->getNumElements();
    _Lglob.resize(nElm);
    _dLdxglob.resize(nElm);
    _dLdyglob.resize(nElm);
    std::vector<double> geoCoord;
    std::vector<double> xq(3, 0.0);

    for(int iElm = 0; iElm < nElm; ++iElm) {
      _mesh->getCoord(_cncGeoTag, iElm, geoCoord);

      std::vector<double> LAtQuadNodes(_nFunctions * _nQuad);
      std::vector<double> dLdxAtQuadNodes(_nFunctions * _nQuad);
      std::vector<double> dLdyAtQuadNodes(_nFunctions * _nQuad);

      for(size_t i = 0; i < _xQuad.size(); ++i) {
        double r[3] = {_xQuad[i], _yQuad[i], _zQuad[i]};
        cnc->getFeSpace()->interpolateVectorField(geoCoord, r, xq);

        std::vector<double> l(_nFunctions, 0.0);
        std::vector<double> dldx(_nFunctions, 0.0);
        std::vector<double> dldy(_nFunctions, 0.0);

        feCheck(Lphys(iElm, xq, l, dldx, dldy));

        for(int j = 0; j < _nFunctions; ++j) LAtQuadNodes[_nFunctions * i + j] = l[j];
        for(int j = 0; j < _nFunctions; ++j) dLdxAtQuadNodes[_nFunctions * i + j] = dldx[j];
        for(int j = 0; j < _nFunctions; ++j) dLdyAtQuadNodes[_nFunctions * i + j] = dldy[j];
        // for(int j = 0; j < _nFunctions; ++j) _dLdt[_nFunctions * i + j] = dldt[j];
      }
      _Lglob[iElm] = LAtQuadNodes;
      _dLdxglob[iElm] = dLdxAtQuadNodes;
      _dLdyglob[iElm] = dLdyAtQuadNodes;
    }
  }

  if(_isGeometricInterpolant) {
    _ignoreNegativeJacobianWarning = false;
    feStatus s = this->getCncGeo()->computeJacobians(_ignoreNegativeJacobianWarning);
    if(s != FE_STATUS_OK) {
      return s;
    }
    this->getCncGeo()->computeConstantElementTransformations();
  }

  return FE_STATUS_OK;
}

/* Computes the inner product of basis functions phi_i, phi_j on element iElm defined by
  (phi_i, phi_j)_iElm = int_iElm phi_i * phi_j ds
  Here we compute all the inner products in a nFunctions x nFunctions Eigen matrix. */
EigenMat feSpace::innerProductBasisFunctions(int iElm)
{
  if(this->getCncGeo()->getDim() == 2) {
    switch(this->getPolynomialDegree()) {
      case 2: {
        EigenMat m = EigenMat::Zero(6, 6);

        const std::vector<double> &J = this->getCncGeo()->getJacobians();

        int n = 6;
        int ex[6] = {0, 1, 0, 2, 1, 0};
        int ey[6] = {0, 0, 1, 0, 1, 2};

        std::cout << _nQuad << std::endl;

        for(int k = 0; k < _nQuad; ++k) {
          for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
              m(i, j) +=
                pow(_xQuad[k], ex[i]) * pow(_yQuad[k], ey[j]) * J[iElm * _nQuad + k] * _wQuad[k];
            }
          }
        }

        return m;
        break;
      }
      default:
        printf("Physical frame basis functions only implemented for polynomials of degree 2.\n");
        exit(-1);
    }
  } else {
    printf("Physical frame basis functions only implemented in 2D.\n");
    exit(-1);
  }
}

double feSpace::innerProductBasisFunctions(int iElm, int ex, int ey)
{
  const std::vector<double> &J = this->getCncGeo()->getJacobians();

  double res = 0.0;
  for(int k = 0; k < _nQuad; ++k) {
    res += pow(_xQuad[k], ex) * pow(_yQuad[k], ey) * J[iElm * _nQuad + k] * _wQuad[k];
  }
  return res;
}

void feSpace::getFunctionsAtQuadNode(const int iQuadNode, std::vector<double> &phi)
{
  for(int i = 0; i < _nFunctions; ++i) {
    phi[i] = _L[_nFunctions * iQuadNode + i];
  }
}

void feSpace::getFunctionsPhysicalGradientAtQuadNode(const int iQuadNode,
                                                     const ElementTransformation &T,
                                                     double *gradPhi)
{
  if(_dim == 0) {
    gradPhi[0] = 0.;
  } else if(_dim == 1) {
    for(int i = 0; i < _nFunctions; ++i) {
      gradPhi[i] = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0];
    }
  } else if(_dim == 2) {
    for(int i = 0; i < _nFunctions; ++i) {
      gradPhi[i * _dim + 0] = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
                              _dLds[_nFunctions * iQuadNode + i] * T.drdx[1];
      gradPhi[i * _dim + 1] = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
                              _dLds[_nFunctions * iQuadNode + i] * T.drdy[1];
    }
  } else if(_dim == 3) {
    for(int i = 0; i < _nFunctions; ++i) {
      gradPhi[i * _dim + 0] = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
                              _dLds[_nFunctions * iQuadNode + i] * T.drdx[1] +
                              _dLdt[_nFunctions * iQuadNode + i] * T.drdx[2];
      gradPhi[i * _dim + 1] = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
                              _dLds[_nFunctions * iQuadNode + i] * T.drdy[1] +
                              _dLdt[_nFunctions * iQuadNode + i] * T.drdy[2];
      gradPhi[i * _dim + 2] = _dLdr[_nFunctions * iQuadNode + i] * T.drdz[0] +
                              _dLds[_nFunctions * iQuadNode + i] * T.drdz[1] +
                              _dLdt[_nFunctions * iQuadNode + i] * T.drdz[2];
    }
  }
}

void feSpace::getVectorFunctionsPhysicalGradientAtQuadNode(const int iQuadNode,
                                                           const ElementTransformation &T,
                                                           std::vector<double> &gradPhi)
{
  assert(gradPhi.size() == (size_t)((unsigned) _dim * _nComponents * _nFunctions));

  int offset = _nFunctions * _nComponents;
  int offsetGrad = _dim * _nComponents;

  if(_dim == 0) {
    // gradPhi[0] = 0.;
    feErrorMsg(FE_STATUS_ERROR, "Not implemented");
    exit(-1);
  } else if(_dim == 1) {
    // for(int i = 0; i < _nFunctions; ++i) {
    //   gradPhi[i] = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0];
    // }
    feErrorMsg(FE_STATUS_ERROR, "Not implemented");
    exit(-1);
  } else if(_dim == 2) {


    for(int i = 0; i < _nFunctions; ++i) {

      // Compute gradient tensor of vector-valued phi
      // Tensor gradient is defined with the convention :
      //      grad(phi)_{mn} = d/dx_m (phi_n)
      double dphi_i_x_dx = _dLdr[offset * iQuadNode + i * _nComponents + 0] * T.drdx[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 0] * T.drdx[1];
      double dphi_i_y_dx = _dLdr[offset * iQuadNode + i * _nComponents + 1] * T.drdx[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 1] * T.drdx[1];
      double dphi_i_x_dy = _dLdr[offset * iQuadNode + i * _nComponents + 0] * T.drdy[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 0] * T.drdy[1];
      double dphi_i_y_dy = _dLdr[offset * iQuadNode + i * _nComponents + 1] * T.drdy[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 1] * T.drdy[1];

      gradPhi[offsetGrad * i + 0] = dphi_i_x_dx;
      gradPhi[offsetGrad * i + 1] = dphi_i_y_dx;
      gradPhi[offsetGrad * i + 2] = dphi_i_x_dy;
      gradPhi[offsetGrad * i + 3] = dphi_i_y_dy;
    }


  } else if(_dim == 3) {
    feErrorMsg(FE_STATUS_ERROR, "Not implemented : 3D gradient of vector shape functions");
    exit(-1);
    // for(int i = 0; i < _nFunctions; ++i) {
    //   gradPhi[i * _dim + 0] = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
    //                           _dLds[_nFunctions * iQuadNode + i] * T.drdx[1] +
    //                           _dLdt[_nFunctions * iQuadNode + i] * T.drdx[2];
    //   gradPhi[i * _dim + 1] = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
    //                           _dLds[_nFunctions * iQuadNode + i] * T.drdy[1] +
    //                           _dLdt[_nFunctions * iQuadNode + i] * T.drdy[2];
    //   gradPhi[i * _dim + 2] = _dLdr[_nFunctions * iQuadNode + i] * T.drdz[0] +
    //                           _dLds[_nFunctions * iQuadNode + i] * T.drdz[1] +
    //                           _dLdt[_nFunctions * iQuadNode + i] * T.drdz[2];
    // }
  }
}

#include "feMatrixInterface.h"

void feSpace::getFunctionsPhysicalHessianAtQuadNode(const int iQuadNode,
                                                    const ElementTransformation &T, double *hessPhi)
{
  // Size hessPhi = nFunctions * [1, 3 or 6]
  if(_dim == 0) {
    hessPhi[0] = 0.;
  } else if(_dim == 1) {
    for(int i = 0; i < _nFunctions; ++i) {
      hessPhi[i] = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[0];
    }
  } else if(_dim == 2) {
    for(int i = 0; i < _nFunctions; ++i) {
      hessPhi[i * 3 + 0] = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[0] +
                           2.0 * _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdx[1];

      hessPhi[i * 3 + 1] = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdy[0] +
                           _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdy[0] +
                           _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdy[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdy[1];

      hessPhi[i * 3 + 2] = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdy[0] * T.drdy[0] +
                           2.0 * _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdy[0] * T.drdy[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdy[1] * T.drdy[1];


      // 2nd order derivatives can also be obtained by solving the
      // following system: (only for linear triangles, otherwise there
      // are additional terms in the RHS, see e.g. 
      // https://scicomp.stackexchange.com/questions/25196/implementing-higher-order-derivatives-for-finite-element
      
      /*
      SquareMatrix mat(3);
      mat(0,0) = pow(T.dxdr[0], 2);
      mat(0,1) = pow(T.dxdr[1], 2);
      mat(0,2) = 2. * T.dxdr[0] * T.dxdr[1];

      mat(1,0) = pow(T.dxds[0], 2);
      mat(1,1) = pow(T.dxds[1], 2);
      mat(1,2) = 2. * T.dxds[0] * T.dxds[1];

      mat(2,0) = T.dxdr[0] * T.dxds[0];
      mat(2,1) = T.dxdr[1] * T.dxds[1];
      mat(2,2) = T.dxdr[0] * T.dxds[1] + T.dxds[0] * T.dxdr[1];

      Vector vec(3);
      vec(0) = _d2Ldr2[_nFunctions * iQuadNode + i];
      vec(1) = _d2Lds2[_nFunctions * iQuadNode + i];
      vec(2) = _d2Ldrs[_nFunctions * iQuadNode + i];

      Vector sol = mat.inverse() * vec;
      feInfo("Comparons:");
      sol.print();
      */
    }
  } else if(_dim == 3) {
    for(int i = 0; i < _nFunctions; ++i) {
      // TODO, one day
      hessPhi[i * 6 + 0] = 0.;
      hessPhi[i * 6 + 1] = 0.;
      hessPhi[i * 6 + 2] = 0.;
      hessPhi[i * 6 + 3] = 0.;
      hessPhi[i * 6 + 4] = 0.;
      hessPhi[i * 6 + 5] = 0.;
    }
  }
}

void feSpace::interpolateField(double *field, int fieldSize, double *r, double *shape, double &res)
{
  res = 0.0;
  L(r, shape);
  feInfo("shape = %f - %f - %f", shape[0], shape[1], shape[2]);
  for(int i = 0; i < fieldSize; ++i) {
    res += field[i] * shape[i];
  }
}

thread_local double SHAPE_FUNCTIONS[18];

double feSpace::interpolateField(std::vector<double> &field, double *r)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    feInfo("Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  L(r, SHAPE_FUNCTIONS);
  // for(int i = 0; i < _nFunctions; ++i) res += field[i] * L(r)[i];
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * SHAPE_FUNCTIONS[i];
  return res;
}

double feSpace::interpolateField(std::vector<double> &field, int iElm, std::vector<double> &x)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * l[i];
  return res;
}

double feSpace::interpolateField(feSolution *sol, std::vector<double> &x)
{
  double u[3];
  int elm = -1;
#ifdef FENG_DEBUG
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x.data(), elm, u);
  if(!isFound) {
    printf(
      "In feSpace::interpolateField : Warning - Point (%f, %f, %f) was not found in the mesh.\n",
      x[0], x[1], x[2]);
    return 0.0;
  }
#else
  static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x.data(), elm, u);
#endif
  std::vector<feInt> adr(this->getNumFunctions());
  std::vector<double> solution(adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  this->initializeAddressingVector(elm, adr);
  for(size_t i = 0; i < adr.size(); ++i) {
    solution[i] = solVec[adr[i]];
  }
  return this->interpolateField(solution, u); // Attention : fonctions de forme locales
}

double feSpace::interpolateField_rDerivative(std::vector<double> &field, double *r) // match
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * dLdr(r)[i];
  return res;
}

double feSpace::interpolateField_sDerivative(std::vector<double> &field, double *r) // match
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * dLds(r)[i];
  return res;
}

double feSpace::interpolateField_xDerivative(std::vector<double> &field, int iElm,
                                             std::vector<double> &x)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * dldx[i];
  return res;
}

double feSpace::interpolateField_yDerivative(std::vector<double> &field, int iElm,
                                             std::vector<double> &x)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * dldy[i];
  return res;
}

void feSpace::interpolateField_gradrs(feSolution *sol, std::vector<double> &x,
                                      std::vector<double> &grad)
{
  double u[3];
  int elm = -1;
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x.data(), elm, u);
  if(!isFound) {
    printf("In feSpace::interpolateField_gradrs: Warning - Point (%f, %f, %f) was not found in the "
           "mesh.\n",
           x[0], x[1], x[2]);
    return;
  } else {
    std::vector<feInt> adr(this->getNumFunctions());
    std::vector<double> solution(adr.size());
    std::vector<double> &solVec = sol->getSolutionReference();
    this->initializeAddressingVector(elm, adr);
    for(size_t i = 0; i < adr.size(); ++i) {
      solution[i] = solVec[adr[i]];
    }
    grad[0] =
      this->interpolateField_rDerivative(solution, u); // Attention : fonctions de forme locales
    grad[1] =
      this->interpolateField_sDerivative(solution, u); // Attention : fonctions de forme locales
  }
}

double feSpace::interpolateFieldAtQuadNode(std::vector<double> &field, int iNode) // match
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _L[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode(std::vector<double> &field, int iElm, int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _Lglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdr[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLds[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_rrDerivative(std::vector<double> &field, int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _d2Ldr2[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_ssDerivative(std::vector<double> &field, int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _d2Lds2[_nFunctions * iNode + i];
  return res;
}

void feSpace::interpolateFieldAtQuadNode_physicalGradient(std::vector<double> &field,
                                                          const int iQuadNode,
                                                          const ElementTransformation &T,
                                                          double *grad)
{
  if(_dim == 0)
  {
    grad[0] = 0.;
  } 
  else if(_dim == 1)
  {
    grad[0] = 0.;
    for(int i = 0; i < _nFunctions; ++i) {
      double dphi_i_dx = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0];
      grad[0] += field[i] * dphi_i_dx;
    }
    // grad[0] /= jac;
  }
  else if(_dim == 2)
  {
    grad[0] = grad[1] = 0.;
    for(int i = 0; i < _nFunctions; ++i) {
      double dphi_i_dx = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
                         _dLds[_nFunctions * iQuadNode + i] * T.drdx[1];
      double dphi_i_dy = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
                         _dLds[_nFunctions * iQuadNode + i] * T.drdy[1];
      grad[0] += field[i] * dphi_i_dx;
      grad[1] += field[i] * dphi_i_dy;
    }
    // grad[0] /= jac;
    // grad[1] /= jac;
  } else if(_dim == 3) {
    grad[0] = grad[1] = grad[2] = 0.;
    for(int i = 0; i < _nFunctions; ++i) {
      double dphi_i_dx = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
                         _dLds[_nFunctions * iQuadNode + i] * T.drdx[1] +
                         _dLdt[_nFunctions * iQuadNode + i] * T.drdx[2];
      double dphi_i_dy = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
                         _dLds[_nFunctions * iQuadNode + i] * T.drdy[1] +
                         _dLdt[_nFunctions * iQuadNode + i] * T.drdy[2];
      double dphi_i_dz = _dLdr[_nFunctions * iQuadNode + i] * T.drdz[0] +
                         _dLds[_nFunctions * iQuadNode + i] * T.drdz[1] +
                         _dLdt[_nFunctions * iQuadNode + i] * T.drdz[2];
      grad[0] += field[i] * dphi_i_dx;
      grad[1] += field[i] * dphi_i_dy;
      grad[2] += field[i] * dphi_i_dz;
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_physicalGradient(const std::vector<double> &field,
                                                                const int nComponents,
                                                                const int iQuadNode,
                                                                const ElementTransformation &T,
                                                                double *grad)
{
  // Size of grad is nComponents x nComponents: derivatives are computed
  // assuming the highest space dimension is the number of vector components.
  // The resulting grad array is:
  //  grad = [ grad(u1), ..., grad(u_nComponents) ]
  for(int i = 0; i < nComponents * nComponents; ++i) grad[i] = 0.;

  int offset = _nFunctions * _nComponents;

  if(nComponents == 1) {
    for(int i = 0; i < _nFunctions; ++i) {
      double dphi_i_dx = _dLdr[offset * iQuadNode + i] * T.drdx[0];
      grad[0] += field[i] * dphi_i_dx;
    }
  } else if(nComponents == 2) {
    for(int i = 0; i < _nFunctions; ++i) {
      // double dphi_i_dx = _dLdr[offset * iQuadNode + i] * T.drdx[0] +
      //                    _dLds[offset * iQuadNode + i] * T.drdx[1];
      // double dphi_i_dy = _dLdr[offset * iQuadNode + i] * T.drdy[0] +
      //                    _dLds[offset * iQuadNode + i] * T.drdy[1];

      // grad[(i % nComponents) * nComponents + 0] += field[i] * dphi_i_dx;
      // grad[(i % nComponents) * nComponents + 1] += field[i] * dphi_i_dy;

      double dphi_i_x_dx = _dLdr[offset * iQuadNode + i * _nComponents + 0] * T.drdx[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 0] * T.drdx[1];
      double dphi_i_y_dx = _dLdr[offset * iQuadNode + i * _nComponents + 1] * T.drdx[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 1] * T.drdx[1];
      double dphi_i_x_dy = _dLdr[offset * iQuadNode + i * _nComponents + 0] * T.drdy[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 0] * T.drdy[1];
      double dphi_i_y_dy = _dLdr[offset * iQuadNode + i * _nComponents + 1] * T.drdy[0]
                         + _dLds[offset * iQuadNode + i * _nComponents + 1] * T.drdy[1];

      grad[0] += field[i] * dphi_i_x_dx;
      grad[1] += field[i] * dphi_i_y_dx;
      grad[2] += field[i] * dphi_i_x_dy;
      grad[3] += field[i] * dphi_i_y_dy;

    }
  } else if(nComponents == 3) {
    for(int i = 0; i < _nFunctions; ++i) {
      feErrorMsg(FE_STATUS_ERROR, "Not implemented");
      exit(-1);
      // double dphi_i_dx = _dLdr[_nFunctions * iQuadNode + i] * T.drdx[0] +
      //                    _dLds[_nFunctions * iQuadNode + i] * T.drdx[1] +
      //                    _dLdt[_nFunctions * iQuadNode + i] * T.drdx[2];

      // double dphi_i_dy = _dLdr[_nFunctions * iQuadNode + i] * T.drdy[0] +
      //                    _dLds[_nFunctions * iQuadNode + i] * T.drdy[1] +
      //                    _dLdt[_nFunctions * iQuadNode + i] * T.drdy[2];

      // double dphi_i_dz = _dLdr[_nFunctions * iQuadNode + i] * T.drdz[0] +
      //                    _dLds[_nFunctions * iQuadNode + i] * T.drdz[1] +
      //                    _dLdt[_nFunctions * iQuadNode + i] * T.drdz[2];

      // FIXME
      // grad[0] += field[i] * dphi_i_dx;
      // grad[1] += field[i] * dphi_i_dy;
      // grad[2] += field[i] * dphi_i_dz;
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_physicalHessian(const std::vector<double> &field,
                                                               const int nComponents,
                                                               const int iQuadNode,
                                                               const ElementTransformation &T,
                                                               double *hessian)
{
  // Size of hessian is nComponents x [1, 3 or 6]: derivatives are computed
  // assuming the highest space dimension is the number of vector components.
  // The resulting hessian array is:
  //  hessian = [ hess(u1), ..., hess(u_nComponents) ]
  // with e.g.
  //  hess(u1) = [u1_xx, u1_xy, u1_yy] in 2D (only unique mixed derivatives)
  int numUniqueMixedDerivatives[3] = {1, 3, 6};
  for(int i = 0; i < numUniqueMixedDerivatives[nComponents-1] * nComponents; ++i)
    hessian[i] = 0.;

  if(nComponents == 1) {
    for(int i = 0; i < _nFunctions; ++i) {
      double d2phi_i_dx2 = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[0];
      hessian[0] += field[i] * d2phi_i_dx2;
    }
  } else if(nComponents == 2) {
    for(int i = 0; i < _nFunctions; ++i) {
      double d2phi_i_dx2 = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[0] +
                           2.0 * _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdx[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdx[1];

      double d2phi_i_dxy = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdy[0] +
                           _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdy[0] +
                           _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdx[0] * T.drdy[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdx[1] * T.drdy[1];

      double d2phi_i_dy2 = _d2Ldr2[_nFunctions * iQuadNode + i] * T.drdy[0] * T.drdy[0] +
                           2.0 * _d2Ldrs[_nFunctions * iQuadNode + i] * T.drdy[0] * T.drdy[1] +
                           _d2Lds2[_nFunctions * iQuadNode + i] * T.drdy[1] * T.drdy[1];

      hessian[(i % nComponents) * 3 + 0] += field[i] * d2phi_i_dx2;
      hessian[(i % nComponents) * 3 + 1] += field[i] * d2phi_i_dxy;
      hessian[(i % nComponents) * 3 + 2] += field[i] * d2phi_i_dy2;
    }
  } else if(nComponents == 3) {
    for(int i = 0; i < _nFunctions; ++i) {
      // TODO, one day
    }
  }
}

double feSpace::interpolateFieldAtQuadNode_xDerivative(std::vector<double> &field, int iElm,
                                                       int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdxglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_yDerivative(std::vector<double> &field, int iElm,
                                                       int iNode)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%ld) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdyglob[iElm][_nFunctions * iNode + i];
  return res;
}

void feSpace::interpolateVectorField(std::vector<double> &field, double *r,
                                     std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
#endif
  std::vector<double> l = L(r);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * l[j];
    }
  }
}

void feSpace::interpolateVectorField(std::vector<double> &field, int nComponents, double *r,
                                     std::vector<double> &res)
{
  // FIXME: rework this, remove nComponents
  // Make this a virtual function of feVectorSpace
  UNUSED(nComponents);
  for(int j = 0; j < _nComponents; ++j) res[j] = 0.0;
  std::vector<double> l = L(r);
  for(int i = 0; i < _nFunctions; ++i) {
    for(int j = 0; j < _nComponents; ++j) {
      res[j] += field[i] * l[i*_nComponents+j];
    }
  }
}

void feSpace::interpolateVectorFieldRT(ElementTransformation &T,
  const int sign,
  const int whichDOF,
  std::vector<double> &field,
  double *r,
  std::vector<double> &res)
{
  UNUSED(whichDOF);
  // for(auto val : field)
  //   feInfo("field = %f", val);
  feInfo("r = [%+-f, %+-f, %+-f]", r[0], r[1], r[2]);
  std::fill(res.begin(), res.end(), 0.);
  std::vector<double> l = L(r);
  // l = {1., 0., 1., 0., 1., 0.};

  double transformed_phi[3] = {0., 0., 0.}, phi[3] = {0., 0., 0.};

  for(int i = 0; i < _nFunctions; ++i) {

    if(i == whichDOF) {

      // Apply contravariant Piola transformation
      for(int j = 0; j < _nComponents; ++j) {
        phi[j] = l[i*_nComponents+j];
      }

      transformed_phi[0] = (T.dxdr[0]*phi[0] + T.dxds[0]*phi[1]) / T.jac;
      transformed_phi[1] = (T.dxdr[1]*phi[0] + T.dxds[1]*phi[1]) / T.jac;

      feInfo("phi  = %+-f - %+-f", phi[0], phi[1]);
      feInfo("tphi = %+-f - %+-f", transformed_phi[0], transformed_phi[1]);

      for(int j = 0; j < _nComponents; ++j) {
        res[j] += sign * field[i] * transformed_phi[j];
      }
    }
  }
}

// Same as interpolateVectorField if _nComponents = 1 and iComponent = 0
// Here the field has the same size as the number of element DOFs, 
// so it contains all the components and we select which one to interpolate.
// This is not the case in interpolateVectorFieldComponentAtQuadnode for now, 
// where the field there has length _nFunctions / _nComponents (only one component)
double feSpace::interpolateVectorFieldComponent(std::vector<double> &field, int iComponent, double *r)
{
  double res = 0.;
  std::vector<double> l = L(r);
  // for(int i = 0; i < _nFunctions / _nComponents; ++i) {
  //   res += field[i * _nComponents + iComponent] * l[i * _nComponents + iComponent];

  //   // _nFunctions * _nComponents * iNode + nComponents * j + i
  // }

  size_t n = field.size();
  for(size_t i = 0; i < n; ++i)
  {
    res += field[i] * l[i * _nComponents + iComponent];
  }
  return res;
}

void feSpace::interpolateVectorField_rDerivative(std::vector<double> &field, double *r,
                                                 std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
#endif
  std::vector<double> dldr = dLdr(r);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * dldr[j];
    }
  }
}

void feSpace::interpolateVectorField_sDerivative(std::vector<double> &field, double *r,
                                                 std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
#endif
  std::vector<double> dlds = dLds(r);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * dlds[j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode(std::vector<double> &field, int iNode,
                                               std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  assert(field.size() == 3 * (unsigned)_nFunctions);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_L[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _L[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode(std::vector<double> &field, int iNode,
                                               std::vector<double> &res, int nComponents)
{
  for(int i = 0; i < nComponents; ++i) res[i] = 0.0;

  // for(int i = 0; i < _nFunctions; ++i) {
  //   res[i % nComponents] += field[i] * _L[_nFunctions * iNode + i];
  // }

  for(int j = 0; j < _nFunctions; ++j) {
    for(int i = 0; i < nComponents; ++i) {
      res[i] += field[j] * _L[_nFunctions * _nComponents * iNode + nComponents * j + i];
    }
  }
}

double feSpace::interpolateVectorFieldComponentAtQuadNode(std::vector<double> &field, int iNode, int iComponent)
{
  UNUSED(iComponent); // (?? Check this at some point)
  feWarning("Component index is not used in interpolation function!");
  double res = 0.;
  for(int i = 0; i < _nFunctions / _nComponents; ++i) {
    res += field[i] * _L[_nFunctions * iNode + _nComponents * i];
  }
  return res;
}

// Same as above but the "field" vector has the same size as the number of DOFs of the space
double feSpace::interpolateVectorFieldComponentAtQuadNode_fullField(std::vector<double> &field, int iNode, int iComponent)
{
  double res = 0.;
  // for(int i = 0; i < _nFunctions / _nComponents; ++i) {
  //   res += field[i * _nComponents + iComponent] * _L[_nFunctions * iNode + _nComponents * i];
  // }
  // return res;

  for(int i = 0; i < _nFunctions; ++i)
  {
    res += field[i] * _L[iNode * _nFunctions * _nComponents +  i * _nComponents + iComponent];
  }
  return res;
}

// Field structure : [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
void feSpace::interpolateVectorFieldAtQuadNode_rDerivative(const std::vector<double> &field, int iNode,
                                                           std::vector<double> &res)
{
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_rDerivative : Erreur - Nombre de valeurs "
           "nodales (%ld) non compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return;
  }
#endif
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * _dLdr[_nFunctions * iNode + j];
    }
  }
}

// Field structure : [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
void feSpace::interpolateVectorFieldAtQuadNode_sDerivative(const std::vector<double> &field, int iNode,
                                                           std::vector<double> &res)
{
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_sDerivative : Erreur - Nombre de valeurs "
           "nodales (%ld) non compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return;
  }
#endif
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _dLds[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_tDerivative(const std::vector<double> &field, int iNode,
                                                           std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    feErrorMsg(FE_STATUS_ERROR,
               "Nombre de valeurs nodales (%ld) non compatible "
               "avec le nombre d'interpolants de l'espace (%d).\n",
               field.size(), (unsigned)_nFunctions);
    exit(-1);
  }
#endif
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * _dLdt[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_rDerivative(const std::vector<double> &field, int iNode,
                                                           double res[3])
{
  res[0] = res[1] = res[2] = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * _dLdr[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_sDerivative(const std::vector<double> &field, int iNode,
                                                           double res[3])
{
  res[0] = res[1] = res[2] = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _dLds[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_tDerivative(const std::vector<double> &field, int iNode,
                                                           double res[3])
{
  res[0] = res[1] = res[2] = 0.0;
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * _dLdt[_nFunctions * iNode + j];
    }
  }
}