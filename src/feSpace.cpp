#include "feSpace.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"

extern int FE_VERBOSE;

feStatus createFiniteElementSpace(feSpace *&space,
  feMesh *mesh,
  const int dimension,
  const geometryType geometry,
  const shapeType shape,
  const int degree,
  const std::string fieldName,
  const std::string cncGeoID,
  const int degreeQuadrature,
  feFunction *fct,
  const bool useGlobalShapeFunctions)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "FINITE ELEMENT SPACE:");
  feInfoCond(FE_VERBOSE > 0, "\t\tCreating FE space for field \"%s\" on entity \"%s\"",
    fieldName.data(),
    cncGeoID.data());
  feInfoCond(FE_VERBOSE > 0, "\t\t\tDimension: %d", dimension);

  if(mesh == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null mesh pointer.");
  if(fct == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null function pointer.");

  if(mesh->getCncGeoTag(cncGeoID) == -1){
    return feErrorMsg(FE_STATUS_ERROR, "Geometric connectivity \"%s\" does not exist in the mesh.",
               cncGeoID.c_str());
  }

  switch(dimension) {

    case 0:
      if(geometry != POINT)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");
      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Point");
      feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Lagrange polynomials");
      space = new feSpace1DP0(mesh, fieldName, cncGeoID, fct);
      break;

    case 1:
      if(geometry != LINE)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Line");

      if(shape == LAGRANGE) {

        feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Lagrange polynomials");
        switch(degree) {
          case 0:
            space = new feSpace1DP0(mesh, fieldName, cncGeoID, fct);
            break;
          case 1:
            space = new feSpace1DP1(mesh, fieldName, cncGeoID, fct);
            break;
          case 2:
            space = new feSpace1DP2(mesh, fieldName, cncGeoID, fct);
            break;
          case 3:
            space = new feSpace1DP3(mesh, fieldName, cncGeoID, fct);
            break;
          case 4:
            space = new feSpace1DP4(mesh, fieldName, cncGeoID, fct);
            break;
          default:
          return feErrorMsg(FE_STATUS_ERROR, 
            "No LAGRANGE 1D finite element space implemented for deg > 4.");
        }

      } else if(shape == LEGENDRE) {

        feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Legendre polynomials");
        space = new feSpace1D_Legendre(degree, mesh, fieldName, cncGeoID, fct);

      } else if(shape == NONCONFORMAL) {

        switch(degree) {
          case 1:
            feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Non-conformal polynomials");
            space = new feSpace1DP1_nonConsistant(mesh, fieldName, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
              "No NONCONFORMAL 1D finite element space implemented for deg > 1.");
        }

      } else {

        return feErrorMsg(FE_STATUS_ERROR, "Unsupported shape functions.");

      }
      break;

    case 2:
      if(geometry != TRI)
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");

      feInfoCond(FE_VERBOSE > 0, "\t\t\tGeometry: Triangle");

      if(shape == LAGRANGE) {
  
        feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Lagrange polynomials");

        switch(degree) {
          case 1:
            space = new feSpaceTriP1(mesh, fieldName, cncGeoID, fct);
            break;
          case 2:
            space = new feSpaceTriP2(mesh, fieldName, cncGeoID, fct, useGlobalShapeFunctions);
            break;
          case 3:
            space = new feSpaceTriP3(mesh, fieldName, cncGeoID, fct);
            break;
          case 4:
            space = new feSpaceTriP4(mesh, fieldName, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
              "No LAGRANGE 2D finite element space implemented on triangles for deg > 4 or = 0.");
        }

      } else if(shape == NONCONFORMAL) {

        feInfoCond(FE_VERBOSE > 0, "\t\t\tShape functions: Non-conformal Crouzeix-Raviart polynomials");

        switch(degree) {
          case 1:
            space = new feSpaceTriP1_nonConsistant(mesh, fieldName, cncGeoID, fct);
            break;
          case 2:
            space = new feSpaceTriP2_nonConsistant(mesh, fieldName, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No TRI_CR finite element space implemented for deg > 2 or = 0.");
        }
      }
      break;

    case 3:
      return feErrorMsg(FE_STATUS_ERROR, "Finite element space in dim = 3 not implemented yet.");

    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Cannot create a finite element space for dimension > 3 or < 0.");
  }

  feInfoCond(FE_VERBOSE > 0, "\t\t\tPolynomial degree: %d", degree);

  // space->useGlobalFunctions(useGlobalShapeFunctions);

  // Set the quadrature rule on this space and the corresponding geometric interpolation space
  feQuadrature rule(degreeQuadrature, dimension, mesh->getCncGeoByName(cncGeoID)->getForme());
  
  // FIXME: the FE space of the geometric connectivity is the one of the last created space...
  feCheck(space->getCncGeo()->getFeSpace()->setQuadratureRule(&rule));
  feCheck(space->setQuadratureRule(&rule));

  return FE_STATUS_OK;
}

feSpace::feSpace(feMesh *mesh, const std::string &fieldID, const std::string &cncGeoID, 
  feFunction *fct, bool useGlobalShapeFunctions)
  : _mesh(mesh)
  , _fieldID(fieldID)
  , _fieldTag(-1)
  , _cncGeoID(cncGeoID)
  , _cncGeoTag(-1)
  , _nQuad(-1)
  , _nFunctions(0)
  , _fct(fct)
  , _useGlobalShapeFunctions(useGlobalShapeFunctions)
{
  // Called with mesh = nullptr when creating the geometric interpolation space
  if(mesh != nullptr){
    _cncGeoTag = mesh->getCncGeoTag(cncGeoID);
  }
}

feCncGeo *feSpace::getCncGeo() {return _mesh->getCncGeoByTag(_cncGeoTag); }

int feSpace::getDim() { return this->getCncGeo()->getDim(); }
int feSpace::getNbElm() { return _mesh->getNbElm(_cncGeoTag); }
int feSpace::getNbNodePerElem() { return _mesh->getNbNodePerElem(_cncGeoTag); }

feStatus feSpace::setQuadratureRule(feQuadrature *rule)
{
  _nQuad = rule->getNQuad();
  _wQuad = rule->getWeights();
  _xQuad = rule->getXPoints();
  _yQuad = rule->getYPoints();
  _zQuad = rule->getZPoints();

  _L.resize(_nFunctions * _nQuad, 0.0);
  _dLdr.resize(_nFunctions * _nQuad, 0.0);
  _dLds.resize(_nFunctions * _nQuad, 0.0);
  _dLdt.resize(_nFunctions * _nQuad, 0.0);

  /* Reference frame discretization : shape functions are computed once on the reference element,
  then evaluated at the quadrature nodes. */
  for(size_t i = 0; i < _xQuad.size(); ++i) {
    double r[3] = {_xQuad[i], _yQuad[i], _zQuad[i]};
    std::vector<double> l = L(r);
    std::vector<double> dldr = dLdr(r);
    std::vector<double> dlds = dLds(r);
    std::vector<double> dldt = dLdt(r);
    for(int j = 0; j < _nFunctions; ++j) _L[_nFunctions * i + j] = l[j];
    for(int j = 0; j < _nFunctions; ++j) _dLdr[_nFunctions * i + j] = dldr[j];
    for(int j = 0; j < _nFunctions; ++j) _dLds[_nFunctions * i + j] = dlds[j];
    for(int j = 0; j < _nFunctions; ++j) _dLdt[_nFunctions * i + j] = dldt[j];
  }

  if(_fieldID != "GEO" && _useGlobalShapeFunctions) {
    feInfo("USING GLOBAL FUNCTIONS on space %s - %s", this->_fieldID.c_str(),
           this->_cncGeoID.c_str());
    feInfo("USING GLOBAL FUNCTIONS");
    feInfo("USING GLOBAL FUNCTIONS");
    /* Physical frame discretization : shape functions are computed on the physical element. */
    feCncGeo *cnc = this->getCncGeo();
    int nElm = cnc->getNbElm();
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

  // If the space is a geometric interpolant, precompute jacobians of the elements
  if(_fieldID == "GEO") {
    feStatus s = this->getCncGeo()->computeJacobians();
    if(s != FE_STATUS_OK){
      return s;
    }
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

        std::vector<double> &J = this->getCncGeo()->getJacobians();

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
  std::vector<double> &J = this->getCncGeo()->getJacobians();

  double res = 0.0;
  for(int k = 0; k < _nQuad; ++k) {
    res += pow(_xQuad[k], ex) * pow(_yQuad[k], ey) * J[iElm * _nQuad + k] * _wQuad[k];
  }
  return res;
}

void feSpace::interpolateField(double *field, int fieldSize, double *r, double *shape, double &res)
{
  res = 0.0;
  L(r, shape);
  for(int i = 0; i < fieldSize; ++i){ res += field[i] * shape[i]; }
}

double feSpace::interpolateField(std::vector<double> &field, double *r)
{
  double res = 0.0;
#ifdef FENG_DEBUG
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
#endif
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * L(r)[i];
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

double feSpace::interpolateField(feNumber *number, feSolution *sol, std::vector<double> &x)
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
  std::vector<feInt> adr(this->getNbFunctions());
  std::vector<double> solution(adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  this->initializeAddressingVector(number, elm, adr);
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

void feSpace::interpolateField_gradrs(feNumber *number, feSolution *sol, std::vector<double> &x,
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
    std::vector<feInt> adr(this->getNbFunctions());
    std::vector<double> solution(adr.size());
    std::vector<double> &solVec = sol->getSolutionReference();
    this->initializeAddressingVector(number, elm, adr);
    for(size_t i = 0; i < adr.size(); ++i) {
      solution[i] = solVec[adr[i]];
    }
    grad[0] = this->interpolateField_rDerivative(
      solution, u); // Attention : fonctions de forme locales
    grad[1] = this->interpolateField_sDerivative(
      solution, u); // Attention : fonctions de forme locales
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

double feSpace::interpolateFieldAtQuadNode_rDerivative(std::vector<double> &field,
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
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdr[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_sDerivative(std::vector<double> &field,
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
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLds[_nFunctions * iNode + i];
  return res;
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
  std::vector<double> l=L(r);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*L(r)[j];
      res[i] += field[3 * j + i] * l[j];
    }
  }
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
  std::vector<double> dldr=dLdr(r);
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
  std::vector<double> dlds=dLds(r);
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
#ifdef FENG_DEBUG
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
#endif
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_L[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _L[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> &field, int iNode,
                                                           std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
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
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _dLdr[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> &field, int iNode,
                                                           std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
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