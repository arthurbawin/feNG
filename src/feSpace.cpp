#include "feSpace.h"
#include "feSpaceTriangle.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feSolution.h"
#include "omp.h"

feStatus createFiniteElementSpace(feSpace *&space, feMesh *mesh, int dim, elemType type, int deg,
                                  std::string fieldID, std::string cncGeoID, int dQuad,
                                  feFunction *fct, bool useGlobalShapeFunctions)
{
  if(mesh == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null mesh pointer.");
  if(fct == nullptr) return feErrorMsg(FE_STATUS_ERROR, "Null function pointer.");

  switch(dim) {
    case 0:
      if(type == POINT) {
        space = new feSpace1DP0(mesh, fieldID, cncGeoID, fct);
      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");
      }
      break;
    case 1:
      if(type == LINE) {
        switch(deg) {
          case 0:
            space = new feSpace1DP0(mesh, fieldID, cncGeoID, fct);
            break;
          case 1:
            space = new feSpace1DP1(mesh, fieldID, cncGeoID, fct);
            break;
          case 2:
            space = new feSpace1DP2(mesh, fieldID, cncGeoID, fct);
            break;
          case 3:
            space = new feSpace1DP3(mesh, fieldID, cncGeoID, fct);
            break;
          case 4:
            space = new feSpace1DP4(mesh, fieldID, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No LINE finite element space implemented for deg > 4.");
        }
      } else if(type == LINE_CR) {
        switch(deg) {
          case 1:
            space = new feSpace1DP1_nonConsistant(mesh, fieldID, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No LINE finite element space implemented for deg > 4.");
        }
      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");
      }
      break;
    case 2:
      if(type == TRI) {
        switch(deg) {
          case 0:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No TRI finite element space implemented for deg > 4.");
          case 1:
            space = new feSpaceTriP1(mesh, fieldID, cncGeoID, fct);
            break;
          case 2:
            space = new feSpaceTriP2(mesh, fieldID, cncGeoID, fct, useGlobalShapeFunctions);
            break;
          case 3:
            space = new feSpaceTriP3(mesh, fieldID, cncGeoID, fct);
            break;
          case 4:
            space = new feSpaceTriP4(mesh, fieldID, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No LINE finite element space implemented for deg > 4.");
        }
      } else if(type == TRI_CR) {
        switch(deg) {
          case 0:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No TRI_CR finite element space implemented for deg >2.");
          case 1:
            space = new feSpaceTriP1_nonConsistant(mesh, fieldID, cncGeoID, fct);
            break;
          case 2:
            space = new feSpaceTriP2_nonConsistant(mesh, fieldID, cncGeoID, fct);
            break;
          default:
            return feErrorMsg(FE_STATUS_ERROR,
                              "No LINE finite element space implemented for deg > 4.");
        }
      } else {
        return feErrorMsg(FE_STATUS_ERROR, "Unsupported geometry.");
      }
      break;

    case 3:
      return feErrorMsg(FE_STATUS_ERROR, "Finite element space in dim = 3 not implemented yet.");
    default:
      return feErrorMsg(FE_STATUS_ERROR,
                        "Cannot create a finite element space for dimension > 3 or < 0.");
  }

  // space->useGlobalFunctions(useGlobalShapeFunctions);

  // Set the quadrature rule on this space and the corresponding geometric interpolation space
  feQuadrature rule(dQuad, dim, mesh->getCncGeoByName(cncGeoID)->getForme());
  feCheck(space->getCncGeo()->getFeSpace()->setQuadratureRule(&rule));
  feCheck(space->setQuadratureRule(&rule));

  return FE_STATUS_OK;
}

feSpace::feSpace(feMesh *mesh, std::string fieldID, std::string cncGeoID, feFunction *fct,
                 bool useGlobalShapeFunctions)
  : _mesh(mesh), _fieldID(fieldID), _fieldTag(-1), _cncGeoID(cncGeoID), _cncGeoTag(-1), _nQuad(-1),
    _nFunctions(0), _fct(fct), _useGlobalShapeFunctions(useGlobalShapeFunctions)
{
  if(mesh != nullptr) {
    _cncGeoTag = mesh->getCncGeoTag(cncGeoID);
    if(_cncGeoTag == -1) {
      feErrorMsg(FE_STATUS_ERROR, "Geometric connectivity \"%s\" does not exist in the mesh.",
                 cncGeoID.c_str());
      std::exit(1);
    }
  }
}

feCncGeo *feSpace::getCncGeo() { return _mesh->getCncGeoByTag(_cncGeoTag); }

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
    double rc[3] = {1. / 3., 1. / 3., 1. / 3.};

    for(int iElm = 0; iElm < nElm; ++iElm) {
      geoCoord = _mesh->getCoord(_cncGeoTag, iElm);

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
    this->getCncGeo()->computeJacobians();
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

// Initialize solution vector from a feSolution
void feSpace::initializeSolution(feSolution *sol)
{
  _sol.resize(_adr.size());
  for(size_t i = 0; i < _adr.size(); ++i) {
    _sol[i] = sol->getSolAtDOF(_adr[i]);
  }
}

// Initialize solution vector from another vector (typically coming from a solutionContainer)
void feSpace::initializeSolution(std::vector<double> &sol)
{
  _sol.resize(_adr.size());
  for(size_t i = 0; i < _adr.size(); ++i) {
    _sol[i] = sol[_adr[i]];
  }
}

void feSpace::initializeSolutionDot(feSolution *sol)
{
  _soldot.resize(_adr.size());
  for(size_t i = 0; i < _adr.size(); ++i) _soldot[i] = sol->getSolDotAtDOF(_adr[i]);
}

double feSpace::interpolateSolution(double r[3])
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * L(r)[i];
  return res;
}

double feSpace::interpolateSolution_rDerivative(double r[3])
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * dLdr(r)[i];
  return res;
}

double feSpace::interpolateSolution_sDerivative(double r[3])
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * dLds(r)[i];
  return res;
}

double feSpace::interpolateSolution(int iElm, std::vector<double> &x)
{
  double res = 0.0;
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * l[i];
  return res;
}

double feSpace::interpolateSolution(feNumber *number, feSolution *sol, std::vector<double> &x)
{
  std::vector<double> u(3, 0.0);
  int elm = -1;
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, u);
  if(!isFound) {
    printf(
      "In feSpace::interpolateSolution : Warning - Point (%f, %f, %f) was not found in the mesh.\n",
      x[0], x[1], x[2]);
    return 0.0;
  } else {
    this->initializeAddressingVector(number, elm);
    this->initializeSolution(sol);
    return this->interpolateSolution(u.data()); // Attention : fonctions de forme locales
  }
}

void feSpace::interpolateSolution_gradrs(feNumber *number, feSolution *sol, std::vector<double> &x,
                                         std::vector<double> &grad)
{
  std::vector<double> u(3, 0.0);
  int elm = -1;
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, u);
  if(!isFound) {
    printf(
      "In feSpace::interpolateSolution : Warning - Point (%f, %f, %f) was not found in the mesh.\n",
      x[0], x[1], x[2]);
    return;
  } else {
    this->initializeAddressingVector(number, elm);
    this->initializeSolution(sol);
    grad[0] =
      this->interpolateSolution_rDerivative(u.data()); // Attention : fonctions de forme locales
    grad[0] =
      this->interpolateSolution_sDerivative(u.data()); // Attention : fonctions de forme locales
  }
}

double feSpace::interpolateSolution_xDerivative(int iElm, std::vector<double> &x)
{
  double res = 0.0;
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * dldx[i];
  return res;
}

double feSpace::interpolateSolution_yDerivative(int iElm, std::vector<double> &x)
{
  double res = 0.0;
  std::vector<double> l(_nFunctions, 0.0);
  std::vector<double> dldx(_nFunctions, 0.0);
  std::vector<double> dldy(_nFunctions, 0.0);
  Lphys(iElm, x, l, dldx, dldy);
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * dldy[i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode(int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) {
    res += _sol[i] * _L[_nFunctions * iNode + i];
    // std::cout<<"_sol [ "<< i<<"]="<<_sol[i] << "_L =" << _L[_nFunctions * iNode + i]<<std::endl;
  }
  return res;
}

double feSpace::interpolateSolutionAtQuadNode(int iElm, int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * _Lglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateSolutionDotAtQuadNode(int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _soldot[i] * _L[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_rDerivative(int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * _dLdr[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_sDerivative(int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * _dLds[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_xDerivative(int iElm, int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * _dLdxglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateSolutionAtQuadNode_yDerivative(int iElm, int iNode)
{
  double res = 0.0;
  for(int i = 0; i < _nFunctions; ++i) res += _sol[i] * _dLdyglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateField(std::vector<double> field, double r[3])
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * L(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode(std::vector<double> field, int iNode)
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _L[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateField_rDerivative(std::vector<double> field, double r[3])
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateField : Erreur - Nombre de valeurs nodales non compatible avec "
           "le nombre d'interpolants de l'espace.\n");
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * dLdr(r)[i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode)
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%d) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdr[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_sDerivative(std::vector<double> field, int iNode)
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%d) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLds[_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_xDerivative(std::vector<double> field, int iElm, int iNode)
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%d) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdxglob[iElm][_nFunctions * iNode + i];
  return res;
}

double feSpace::interpolateFieldAtQuadNode_yDerivative(std::vector<double> field, int iElm, int iNode)
{
  double res = 0.0;
  if(field.size() != (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateFieldAtQuadNode : Erreur - Nombre de valeurs nodales (%d) non "
           "compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return res;
  }
  for(int i = 0; i < _nFunctions; ++i) res += field[i] * _dLdyglob[iElm][_nFunctions * iNode + i];
  return res;
}

void feSpace::interpolateVectorField(std::vector<double> field, double r[3],
                                     std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*L(r)[j];
      res[i] += field[3 * j + i] * L(r)[j];
    }
  }
}

void feSpace::interpolateVectorField_rDerivative(std::vector<double> field, double r[3],
                                                 std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * dLdr(r)[j];
    }
  }
}

void feSpace::interpolateVectorField_sDerivative(std::vector<double> field, double r[3],
                                                 std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorField : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      res[i] += field[3 * j + i] * dLds(r)[j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode(std::vector<double> field, int iNode,
                                               std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode : Erreur - Nombre de valeurs nodales non "
           "compatible avec le nombre d'interpolants de l'espace.\n");
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_L[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _L[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_rDerivative(std::vector<double> field, int iNode,
                                                           std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_rDerivative : Erreur - Nombre de valeurs "
           "nodales (%d) non compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _dLdr[_nFunctions * iNode + j];
    }
  }
}

void feSpace::interpolateVectorFieldAtQuadNode_sDerivative(std::vector<double> field, int iNode,
                                                           std::vector<double> &res)
{
  // Field structure :
  // [fx0 fy0 fz0 fx1 fy1 fz1 ... fxn fyn fzn]
  res[0] = res[1] = res[2] = 0.0;
  if(field.size() != 3 * (unsigned)_nFunctions) {
    printf(" In feSpace::interpolateVectorFieldAtQuadNode_sDerivative : Erreur - Nombre de valeurs "
           "nodales (%d) non compatible avec le nombre d'interpolants de l'espace (%d).\n",
           field.size(), (unsigned)_nFunctions);
    return;
  }
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < _nFunctions; ++j) {
      // res[i] += field[3*i+j]*_dLdr[_nFunctions*iNode+j];
      res[i] += field[3 * j + i] * _dLds[_nFunctions * iNode + j];
    }
  }
}

void feSpace::printL()
{
  std::cout << "taille = " << _L.size() << std::endl;
  for(double l : _L) std::cout << l << " ";
  std::cout << std::endl;
}

void feSpace::printdLdr()
{
  std::cout << "taille = " << _dLdr.size() << std::endl;
  for(double l : _dLdr) std::cout << l << " ";
  std::cout << std::endl;
}

void feSpace1DP0::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) number->defDDLSommet(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
}

void feSpace1DP0::initializeAddressingVector(feNumber *number, int numElem)
{
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
}

void feSpace1DP0::initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr)
{
  // for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i)
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
}

void feSpace1DP1::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
  }
}

void feSpace1DP1::initializeAddressingVector(feNumber *number, int numElem)
{
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
}

void feSpace1DP1::initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
}

void feSpace1DP1_nonConsistant::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLElement(_mesh, _cncGeoID, i, 1);
    // If the line is a boundary element, the edge should be set by the interior element
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
  }
}

void feSpace1DP1_nonConsistant::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    // Set essential BC on the edge if the line is a boundary element
    // If the line is an interior element, there is nothing in cncGeo->connecEdges
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP1_nonConsistant::initializeAddressingVector(feNumber *number, int numElem)
{
  _adr[0] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
}

void feSpace1DP1_nonConsistant::initializeAddressingVector(feNumber *number, int numElem,
                                                           std::vector<int> &adr)
{
  adr[0] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
}

void feSpace1DP2::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // Loop over the elements nodes on geometric interpolant (> 2 if curved)
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet(_mesh, _cncGeoID, i, j);
    }
    // If the edge is curved (P2), the middle vertex is already numbered.
    // This is true for Pn geometries with n even, where the middle vertex matches
    // the middle node of the interpolant, but so far we will limit to P2 geometries.
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLElement(_mesh, _cncGeoID, i, 1);
      // If the line is a boundary element, the edge should be set by the interior element

      // FIXME : When should this be commented ? :/
      number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
    }
  }
}

void feSpace1DP2::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    for(int j = 0; j < this->getCncGeo()->getNbNodePerElem(); ++j) {
      number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, j);
    }
    if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() != 2) {
      number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
      // Set essential BC on the edge if the line is a boundary element
      // If the line is an interior element, there is nothing in cncGeo->connecEdges
      // This is only true is there is a surface feSpace associated with the current feSpace !

      number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
    }
  }
}

void feSpace1DP2::initializeAddressingVector(feNumber *number, int numElem)
{
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    _adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  } else {
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

void feSpace1DP2::initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  if(this->getCncGeo()->getFeSpace()->getPolynomialDegree() == 2) {
    adr[2] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

void feSpace1DP3::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // Do not loop over all geometric vertices : only 0 and 1
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    /* The mid-edge vertex does not match any of the element dofs,
    so they must both be added. */
    number->defDDLElement(_mesh, _cncGeoID, i, 2);
    // TODO : add a test to check if it is in the situation of Verwer
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 2);
  }
}

void feSpace1DP3::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP3::initializeAddressingVector(feNumber *number, int numElem)
{
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  } else {
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

void feSpace1DP3::initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
  } else {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
  }
}

void feSpace1DP4::initializeNumberingUnknowns(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    number->defDDLSommet(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet(_mesh, _cncGeoID, i, 1);
    number->defDDLElement(_mesh, _cncGeoID, i, 3);
    // TODO : add a test to check if it is in the situation of Verwer
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 1);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 2);
    number->defDDLEdge(_mesh, _cncGeoID, i, 0, 3);
  }
}

void feSpace1DP4::initializeNumberingEssential(feNumber *number)
{
  for(int i = 0; i < _mesh->getNbElm(_cncGeoID); ++i) {
    // TODO : Modifier pour elements courbes
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 0);
    number->defDDLSommet_essentialBC(_mesh, _cncGeoID, i, 1);
    number->defDDLElement_essentialBC(_mesh, _cncGeoID, i);
    number->defDDLEdge_essentialBC(_mesh, _cncGeoID, i, 0);
  }
}

void feSpace1DP4::initializeAddressingVector(feNumber *number, int numElem)
{
  _adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  _adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    _adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  } else {
    _adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    _adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    _adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  }
}

void feSpace1DP4::initializeAddressingVector(feNumber *number, int numElem, std::vector<int> &adr)
{
  adr[0] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 0);
  adr[1] = number->getDDLSommet(_mesh, _cncGeoID, numElem, 1);
  int e = _mesh->getEdge(_cncGeoID, numElem, 0);
  if(e > 0) {
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  } else {
    adr[4] = number->getDDLElement(_mesh, _cncGeoID, numElem, 0);
    adr[3] = number->getDDLElement(_mesh, _cncGeoID, numElem, 1);
    adr[2] = number->getDDLElement(_mesh, _cncGeoID, numElem, 2);
  }
}