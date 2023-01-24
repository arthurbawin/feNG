#include "feNorm.h"

feStatus createNorm(feNorm *&norm,
  normType type,
  const std::vector<feSpace*> &spaces,
  feSolution *sol,
  feFunction *scalarSolution,
  feVectorFunction *vectorSolution)
{
  // Check that all interpolation spaces are defined on the same connectivity
  int cncGeoTag = spaces[0]->getCncGeoTag();
  for(size_t i = 0; i < spaces.size(); ++i) {
    if(spaces[i]->getCncGeoTag() != cncGeoTag)
      return feErrorMsg(FE_STATUS_ERROR, "Norm defined on more than one geometric connectivity (physical entity).\n"
        "All finite element spaces in \"spaces\" should be defined on the same geometric connectivity.\n");
  }

  norm = new feNorm(type, spaces, sol, scalarSolution, vectorSolution);

  return FE_STATUS_OK;
}

feNorm::feNorm(normType type, const std::vector<feSpace*> &spaces, feSolution *sol, feFunction *scalarSolution,
               feVectorFunction *vectorSolution)
  : _type(type)
  , _spaces(spaces)
  , _solution(sol)
  , _cnc(spaces[0]->getCncGeo())
  , _scalarSolution(scalarSolution)
  , _vectorSolution(vectorSolution)
  , _J(_cnc->getJacobians())
{
  _nQuad = spaces[0]->getNumQuadPoints();
  _w = spaces[0]->getQuadratureWeights();
  _geoSpace = _cnc->getFeSpace();
  _nElm = _spaces[0]->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * _cnc->getNumVerticesPerElem());

  _adr.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k)
    _adr[k].resize(spaces[k]->getNumFunctions());
  _localSol.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k)
    _localSol[k].resize(spaces[k]->getNumFunctions());

  _cncOnly = false;
}

feNorm::feNorm(feCncGeo *cnc, feFunction *scalarSolution, feVectorFunction *vectorSolution)
  : _cnc(cnc)
  , _scalarSolution(scalarSolution)
  , _vectorSolution(vectorSolution)
  , _J(_cnc->getJacobians())
{
  _nQuad = cnc->getFeSpace()->getNumQuadPoints();
  _w = cnc->getFeSpace()->getQuadratureWeights();
  _geoSpace = cnc->getFeSpace();
  _nElm = cnc->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * cnc->getNumVerticesPerElem());

  _cncOnly = true;
}

double feNorm::compute(normType type)
{
  if( _cncOnly && !(type == AREA || type == INTEGRAL_F) )
  {
    feErrorMsg(FE_STATUS_ERROR, "This norm was created with a geometric"
      " connectivity only. Cannot compute this norm without"
      " providing solution and finite element space.");
    exit(-1);
  }

  double res;
  switch(type)
  {
    case L1:
      res = this->computeLpNorm(1, false);
      break;
    case L1_ERROR:
      res = this->computeLpNorm(1, true);
      break;
    case L2:
      res = this->computeLpNorm(2, false);
      break;
    case L2_ERROR:
      res = this->computeLpNorm(2, true);
      break;
    case LINF:
      res = this->computeLInfNorm(false);
      break;
    case LINF_ERROR:
      res = this->computeLInfNorm(true);
      break;
    case AREA:
      res = this->computeArea();
      break;
    case INTEGRAL:
      return this->computeIntegral();
      break;
    case INTEGRAL_F:
      return this->computeIntegralUserFunction();
      break;
    case DOT_PRODUCT:
      return this->computeIntegralDotProduct();
      break;
  }
  return res;
}

double feNorm::compute()
{
  return this->compute(_type);
}

void feNorm::initializeLocalSolutionOnSpace(int iSpace, int iElm)
{
  _spaces[iSpace]->initializeAddressingVector(iElm, _adr[iSpace]);
  _solution->getSolAtDOF(_adr[iSpace], _localSol[iSpace]);
}

double feNorm::computeLpNorm(int p, bool error)
{
  double res = 0.0, uh, u, t = _solution->getCurrentTime();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {

    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = error ? _scalarSolution->eval(t, _pos) : 0.0;

      res += pow(fabs(u - uh), p) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1./(double) p);
}

double feNorm::computeL1Norm(bool error)
{
  double res = 0.0, uh, u, t = _solution->getCurrentTime();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {

    this->initializeLocalSolutionOnSpace(0, iElm);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = error ? _scalarSolution->eval(t, _pos) : 0.0;
      res += fabs(u - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeL2Norm(bool error)
{
  double res = 0.0, uh, u, t = _solution->getCurrentTime();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {

    this->initializeLocalSolutionOnSpace(0, iElm);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = error ? _scalarSolution->eval(t, _pos) : 0.0;
      res += (u - uh) * (u - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return sqrt(res);
}

double feNorm::computeLInfNorm(bool error)
{
  double res = -DBL_MAX, uh, u, t = _solution->getCurrentTime();

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {

    this->initializeLocalSolutionOnSpace(0, iElm);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = error ? _scalarSolution->eval(t, _pos) : 0.0;
      res = fmax(res, fabs(u - uh));
    }
  }
  return res;
}

double feNorm::computeH1Norm(bool error)
{
  double res = 0.0;

  /* ... */
  
  return res;
}

double feNorm::computeArea()
{
  double res = 0.0;
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    for(int k = 0; k < _nQuad; ++k) {
      res += _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeIntegral()
{
  double res = 0.0, uh;
  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    this->initializeLocalSolutionOnSpace(0, iElm);
    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      res += uh * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeIntegralUserFunction()
{
  double res = 0.0, u, t = 0.;

  #if defined(HAVE_OMP)
  #pragma omp parallel for private(_geoCoord, u) reduction(+ : res) schedule(dynamic)
  #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {

    _geoSpace->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = (*_scalarSolution)(t, _pos);
      res += u * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

// Compute the integral of uh * f, where f is a user-defined function 
double feNorm::computeIntegralDotProduct()
{
  // double I = 0.0, solInt, J, t = _solution->getCurrentTime();
  
  // // std::vector<feInt> adr0(_spaces[0]->getNumFunctions());
  // // std::vector<double> sol0(adr0.size());
  // // std::vector<double> &solVec = _solution->getSolutionReference();

  // for(int iElm = 0; iElm < _nElm; ++iElm) {


  //   _spaces[0]->initializeAddressingVector(iElm, adr0);
  //   for(size_t i = 0; i < adr0.size(); ++i) {
  //     sol0[i] = solVec[adr0[i]];
  //   }

  //   _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);
  //   for(int k = 0; k < _nQuad; ++k) {
  //     solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
  //     std::vector<double> x(3, 0.0);
  //     _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
  //     std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
  //     std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
  //     _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k, dxdr);
  //     _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k, dxds);
  //     J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
  //     I += _scalarSolution->eval(t, x) * solInt * J * _w[k];
  //   }
  // }
  // return I;
}

/* Estimates the L2 norm of the error taking an external solution as the reference solution. */
feStatus feNorm::computeErrorNormFromExternalSolution(feMetaNumber *metaNumber, feSolution *sol,
                                                      feMesh *mesh, feMetaNumber *refMN,
                                                      feSolution *refSol, feMesh *refMesh,
                                                      const std::vector<feSpace *> refSpaces,
                                                      double &res)
{
  double normL2 = 0.0, solInt, solRef;
  std::vector<feInt> adr(_spaces[0]->getNumFunctions());
  std::vector<double> sol0(adr.size());
  std::vector<double> &solVec = sol->getSolutionReference();
  std::vector<double> &solVecS = refSol->getSolutionReference();

  if(fabs(sol->getCurrentTime() - refSol->getCurrentTime()) > 1e-3) {
    feWarning("Solutions time differ by  more than 1e-3 : tRef "
              "= %f - tSol = %f.\n",
              refSol->getCurrentTime(), sol->getCurrentTime());
    feWarning("Maybe the solution file does not match the current solution.\n");
  }

  std::vector<double> x(3, 0.0);
  double r[3];

  // Identify the feSpace from among the set of spaces from the reference solution (we should use a
  // map) Aussi : la numérotation doit être recréée pour assigner correctement les ddl aux espaces.
  // Donc soit elle est la même que pour le problème courant, ce qui suppose que les fespaces sont
  // les memes, ce qui a du sens, soit il faut redéfinir un nouveau set de feSpace à partir d'une
  // autre simulation, ce qui semble hasardeux.
  bool matchingSpace = false;
  for(feSpace *fS : refSpaces) {
    std::vector<feInt> adrS(fS->getNumFunctions());
    std::vector<double> solS(adrS.size());
    if(fS->getFieldID() == _spaces[0]->getFieldID() &&
       fS->getCncGeoID() == _spaces[0]->getCncGeoID()) {
      matchingSpace = true;

      // feNumber *n = refMN->getNumbering(fS->getFieldID());

      for(int iElm = 0; iElm < _nElm; ++iElm) {
        _spaces[0]->initializeAddressingVector(iElm, adr);
        for(size_t i = 0; i < adr.size(); ++i) {
          sol0[i] = solVec[adr[i]];
        }

        mesh->getCoord(_cnc, iElm, _geoCoord);
        for(int k = 0; k < _nQuad; ++k) {
          _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, x);
          // Current solution
          solInt = _spaces[0]->interpolateFieldAtQuadNode(sol0, k);
          // Reference solution evaluated at quad node
          int elm = -1;
          bool isFound = static_cast<feMesh2DP1 *>(refMesh)->locateVertex(x.data(), elm, r);
          if(!isFound) {
            feWarning("Point (%f, %f, %f) was not found in the mesh.\n", x[0], x[1], x[2]);
            solRef = solInt; // Points outside the mesh do not contribute
          } else {
            fS->initializeAddressingVector(elm, adrS);
            for(size_t i = 0; i < adr.size(); ++i) {
              solS[i] = solVecS[adrS[i]];
            }
            solRef = fS->interpolateField(solS, r);
            normL2 += (solInt - solRef) * (solInt - solRef) * _J[_nQuad * iElm + k] * _w[k];
          }

          // normL2 += (solInt - solRef) * (solInt - solRef) * J[_nQuad * iElm + k] * _w[k];
        }
      }

      res = sqrt(normL2);

      break;
    }
  }

  if(!matchingSpace)
    return feErrorMsg(FE_STATUS_ERROR, "No finite element space in the reference set "
                                       "matches the target space on this connectivity.\n");

  return FE_STATUS_OK;
}
