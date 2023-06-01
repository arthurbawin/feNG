#include "feNorm.h"

thread_local std::vector<double> MY_POS(3, 0.);
thread_local std::vector<double> MY_GEOCOORD(18, 0.);

feStatus createNorm(feNorm *&norm, normType type, const std::vector<feSpace *> &spaces,
                    feSolution *sol, feFunction *scalarSolution, feVectorFunction *vectorSolution)
{
  for(auto *fS : spaces) {
    if(fS == nullptr) {
      return feErrorMsg(FE_STATUS_ERROR,
                        "Null pointer in vector of FE spaces, maybe you forgot to initialize it.");
    }
  }

  // Check that all interpolation spaces are defined on the same connectivity
  int cncGeoTag = spaces[0]->getCncGeoTag();
  for(auto *fS : spaces) {
    if(fS->getCncGeoTag() != cncGeoTag)
      return feErrorMsg(FE_STATUS_ERROR,
                        "Norm defined on more than one geometric connectivity (physical entity).\n"
                        "All finite element spaces in \"spaces\" should be defined on the same "
                        "geometric connectivity.\n");
  }

  norm = new feNorm(type, spaces, sol, scalarSolution, vectorSolution);

  return FE_STATUS_OK;
}

feNorm::feNorm(normType type, const std::vector<feSpace *> &spaces, feSolution *sol,
               feFunction *scalarSolution, feVectorFunction *vectorSolution)
  : _type(type), _spaces(spaces), _solution(sol), _scalarSolution(scalarSolution),
    _vectorSolution(vectorSolution), _cnc(spaces[0]->getCncGeo()), _J(_cnc->getJacobians())
{
  _nQuad = spaces[0]->getNumQuadPoints();
  _w = spaces[0]->getQuadratureWeights();
  _geoSpace = _cnc->getFeSpace();
  _nElm = _spaces[0]->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * _cnc->getNumVerticesPerElem());

  _adr.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k) _adr[k].resize(spaces[k]->getNumFunctions());
  _localSol.resize(spaces.size());
  for(size_t k = 0; k < spaces.size(); ++k) _localSol[k].resize(spaces[k]->getNumFunctions());

  // Provide only a geometric connectivity, to compute integrals of source terms, etc.
  // but not of the finite element solution.
  _cncOnly = false;
}

feNorm::feNorm(feCncGeo *cnc, feFunction *scalarSolution, feVectorFunction *vectorSolution)
  : _cnc(cnc), _scalarSolution(scalarSolution), _vectorSolution(vectorSolution),
    _J(_cnc->getJacobians())
{
  _nQuad = cnc->getFeSpace()->getNumQuadPoints();
  _w = cnc->getFeSpace()->getQuadratureWeights();
  _geoSpace = cnc->getFeSpace();
  _nElm = cnc->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * cnc->getNumVerticesPerElem());

  // Provide only a geometric connectivity, to compute integrals of source terms, etc.
  // but not of the finite element solution.
  _cncOnly = true;
}

double feNorm::compute(normType type)
{
  if(_cncOnly && !(type == AREA || type == INTEGRAL_F)) {
    feErrorMsg(FE_STATUS_ERROR, "This norm was created with a geometric"
                                " connectivity only. Cannot compute this norm without"
                                " providing solution and finite element space.");
    exit(-1);
  }

  double res;
  switch(type) {
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
    case L2_ERROR_ESTIMATE:
      res = this->computeLpErrorEstimator(2);
      break;
    case VECTOR_L2:
      res = this->computeVectorLpNorm(2, false);
      break;
    case VECTOR_L2_ERROR:
      res = this->computeVectorLpNorm(2, true);
      break;
    case LINF:
      res = this->computeLInfNorm(false);
      break;
    case LINF_ERROR:
      res = this->computeLInfNorm(true);
      break;
    case SEMI_H1:
      res = this->computeH1SemiNorm(2, false);
      break;
    case SEMI_H1_ERROR:
      res = this->computeH1SemiNorm(2, true);
      break;
    case SEMI_H1_ERROR_ESTIMATE:
      res = this->computeH1SemiNormErrorEstimator();
      break;
    case H1:
      res = this->computeH1Norm(false);
      break;
    case H1_ERROR:
      res = this->computeH1Norm(true);
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

double feNorm::compute() { return this->compute(_type); }

void feNorm::initializeLocalSolutionOnSpace(int iSpace, int iElm)
{
  _spaces[iSpace]->initializeAddressingVector(iElm, _adr[iSpace]);
  _solution->getSolAtDOF(_adr[iSpace], _localSol[iSpace]);
}

double feNorm::computeLpNorm(int p, bool error)
{
  double res = 0.0, uh, u, t = _solution->getCurrentTime();

  if(error && _scalarSolution == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp norm of error function"
                                " because exact solution is NULL");
    exit(-1);
  }

  // FILE *myfile = fopen("errorOnElements.pos", "w");
  // fprintf(myfile, "View \"errorOnElements\"{\n");

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    double eLocPowP = 0.;
    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      u = error ? _scalarSolution->eval(t, _pos) : 0.0;

      eLocPowP += pow(fabs(u - uh), p) * _J[_nQuad * iElm + k] * _w[k];
    }

    res += eLocPowP;

    double eLoc = pow(eLocPowP, 1. / (double)p);
    // Plot error
    // _cnc->writeElementToPOS(myfile, _geoCoord, eLoc);
  }

  // fprintf(myfile, "};\n"); fclose(myfile);

  return pow(res, 1. / (double)p);
}

double feNorm::computeLpNormOnElement(int p, bool error, int iElm)
{
  double uh, u, t = _solution->getCurrentTime();

  if(error && _scalarSolution == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp norm of error function"
                                " because exact solution is NULL");
    exit(-1);
  }

  this->initializeLocalSolutionOnSpace(0, iElm);
  _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

  double eLocPowP = 0.;
    // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord) reduction(+ : res) schedule(dynamic)
  // #endif
  for(int k = 0; k < _nQuad; ++k) {
    uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
    _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
    u = error ? _scalarSolution->eval(t, _pos) : 0.0;

    eLocPowP += pow(fabs(u - uh), p) * _J[_nQuad * iElm + k] * _w[k];
  }

  return eLocPowP;
}

double area(const std::vector<double> &triCoord)
{
  // radius of circle circumscribing a triangle
  double dist[3], k = 0.0;
  for(int i = 0; i < 3; i++) {
    double x0 = triCoord[3 * i + 0];
    double y0 = triCoord[3 * i + 1];
    double x1 = triCoord[3 * ((i + 1) % 3) + 0];
    double y1 = triCoord[3 * ((i + 1) % 3) + 1];
    dist[i] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    k += 0.5 * dist[i];
  }
  double const area = std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
  return area;
}

double feNorm::computeSquaredErrorOnElement(int iElm)
{
  double t = _solution->getCurrentTime();

  this->initializeLocalSolutionOnSpace(0, iElm);
  _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

  double eLocPowP = 0.;
  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double u, uh;

    #if defined(HAVE_OMP)
    #pragma omp for reduction(+:eLocPowP)
    #endif
    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, MY_POS);
      u = _scalarSolution->eval(t, MY_POS);
      eLocPowP += (u - uh) * (u - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }


  return eLocPowP;
}

double feNorm::computeLpErrorEstimator(int p)
{
  double res = 0.0, uh, uRec;

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      // _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      // uh = _scalarSolution->eval(t, _pos);
      uRec = _rec->evaluateRecoveryAtQuadNode(PPR::RECOVERY, 0, iElm, k);

      res += pow(fabs(uRec - uh), p) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

double feNorm::computeLpErrorExactVsEstimator(int p)
{
  double res = 0.0, uExact, uRec, t = _solution->getCurrentTime();

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      uExact = _scalarSolution->eval(t, _pos);
      uRec = _rec->evaluateRecoveryAtQuadNode(PPR::RECOVERY, 0, iElm, k);

      res += pow(fabs(uRec - uExact), p) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

thread_local std::vector<double> GRAD_U_EXACT(2, 0.);

double feNorm::computeSemiH1ErrorExactVsEstimator(int p)
{
  double res = 0.0, t = _solution->getCurrentTime();
  std::vector<double> graduExact(2, 0.);

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double normSquared, gradRec[2];

    #if defined(HAVE_OMP)
    #pragma omp for reduction (+:res)
    #endif
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, MY_GEOCOORD);

      for(int k = 0; k < _nQuad; ++k) {
        _geoSpace->interpolateVectorFieldAtQuadNode(MY_GEOCOORD, k, MY_POS);
        _vectorSolution->eval(t, MY_POS, GRAD_U_EXACT);
        gradRec[0] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 0, iElm, k);
        gradRec[1] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 1, iElm, k);

        normSquared = (GRAD_U_EXACT[0] - gradRec[0]) * (GRAD_U_EXACT[0] - gradRec[0])
                    + (GRAD_U_EXACT[1] - gradRec[1]) * (GRAD_U_EXACT[1] - gradRec[1]);

        res += normSquared * _J[_nQuad * iElm + k] * _w[k];
      }
    }
  }

  feInfo("Computed norm of gradient error in %f s", toc());

  return sqrt(res);
}

thread_local std::vector<double> HESS_U_EXACT(4, 0.);

double feNorm::computeErrorHessianExactVsEstimator(int p)
{
  double res = 0.0, t = _solution->getCurrentTime();

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double normSquared, hessRec[4];

    #if defined(HAVE_OMP)
    #pragma omp for reduction(+:res)
    #endif
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, MY_GEOCOORD);

      for(int k = 0; k < _nQuad; ++k) {
        _geoSpace->interpolateVectorFieldAtQuadNode(MY_GEOCOORD, k, MY_POS);
        _vectorSolution->eval(t, MY_POS, HESS_U_EXACT);

        // Averaged derivatives-wise
        hessRec[0] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 2, iElm, k);
        hessRec[1] = (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 3, iElm, k) + 
                      _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 4, iElm, k))/2.;
        hessRec[2] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 5, iElm, k); 

        normSquared =      (HESS_U_EXACT[0] - hessRec[0]) * (HESS_U_EXACT[0] - hessRec[0]) + 
                    + 2. * (HESS_U_EXACT[1] - hessRec[1]) * (HESS_U_EXACT[1] - hessRec[1]) +
                         + (HESS_U_EXACT[3] - hessRec[2]) * (HESS_U_EXACT[3] - hessRec[2]);

        // Component-wise
        // hessRec[0] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 2, iElm, k);
        // hessRec[1] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 3, iElm, k);
        // hessRec[2] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 4, iElm, k); 
        // hessRec[3] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 5, iElm, k);

        // normSquared = 0.;
        // for(int i = 0; i < 4; ++i)
        //   normSquared += (HESS_U_EXACT[i] - hessRec[i]) * (HESS_U_EXACT[i] - hessRec[i]);

        res += normSquared * _J[_nQuad * iElm + k] * _w[k];
      }
    }
  }

  feInfo("Computed norm of hessian error in %f s", toc());

  return sqrt(res);
}

thread_local std::vector<double> D3U_EXACT(8, 0.);

double feNorm::computeErrorThirdDerivativesExactVsEstimator(int p)
{
  double res = 0.0, t = _solution->getCurrentTime();

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double normSquared, d3Rec[4];

    #if defined(HAVE_OMP)
    #pragma omp for reduction(+:res)
    #endif
    for(int iElm = 0; iElm < _nElm; ++iElm) {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, MY_GEOCOORD);

      for(int k = 0; k < _nQuad; ++k) {
        _geoSpace->interpolateVectorFieldAtQuadNode(MY_GEOCOORD, k, MY_POS);
        _vectorSolution->eval(t, MY_POS, D3U_EXACT);

        // Averaged derivatives-wise
        // Be careful for the order of mixed derivatives:
        // it's xxx, [xxy, xyx, xyy, yxx, yxy, yyx], yyy
        // so 7-8-10 and 9-11-12 (first 6 positions are 2 first and 4 second derivatives)
        d3Rec[0] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 6, iElm, k);
        d3Rec[1] = (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 7, iElm, k) +
                    _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 8, iElm, k) +  
                    _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 10, iElm, k))/3.; 
        d3Rec[2] = (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 9, iElm, k) +
                    _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 11, iElm, k) +  
                    _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 12, iElm, k))/3.; 
        d3Rec[3] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 13, iElm, k);

        normSquared =  (D3U_EXACT[0] - d3Rec[0]) * (D3U_EXACT[0] - d3Rec[0])
                + 3. * (D3U_EXACT[1] - d3Rec[1]) * (D3U_EXACT[1] - d3Rec[1])
                + 3. * (D3U_EXACT[3] - d3Rec[2]) * (D3U_EXACT[3] - d3Rec[2])
                +      (D3U_EXACT[7] - d3Rec[3]) * (D3U_EXACT[7] - d3Rec[3]) ;

        // // Component-wise
        // normSquared = 0.;
        // for(int ii = 0; ii < 8; ++ii) {
        //   d3Rec[ii] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 6+ii, iElm, k);
        //   normSquared += (D3U_EXACT[ii] - d3Rec[ii]) * (D3U_EXACT[ii] - d3Rec[ii]);
        // }

        res += normSquared * _J[_nQuad * iElm + k] * _w[k];
      }
    }
  }

  feInfo("Computed norm of d3u error in %f s", toc());

  return sqrt(res);
}


double feNorm::computeVectorLpNorm(int p, bool error)
{
  double res = 0.0, compNorm, t = _solution->getCurrentTime();
  std::vector<double> u(3, 0.);
  std::vector<double> uh(3, 0.);
  int nComponents = _spaces[0]->getNumComponents();

  if(error && _vectorSolution == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute Lp norm of vector error function"
                                " because exact solution is NULL");
    exit(-1);
  }

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord) reduction(+ : res) schedule(dynamic)
  // #endif
  for(int iElm = 0; iElm < _nElm; ++iElm) {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(int k = 0; k < _nQuad; ++k) {
      _spaces[0]->interpolateVectorFieldAtQuadNode(_localSol[0], k, uh, nComponents);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      if(error) (*_vectorSolution)(t, _pos, u);

      compNorm = 0.;
      for(int i = 0; i < nComponents; ++i) {
        compNorm += pow(fabs(u[i] - uh[i]), p);
      }

      res += compNorm * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

double feNorm::computeL1Norm(bool error)
{
  double res = 0.0, uh, u, t = _solution->getCurrentTime();

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  // #endif
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

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  // #endif
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

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res) schedule(dynamic)
  // #endif
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

double feNorm::computeH1SemiNorm(int p, bool error)
{
  double res = 0.0, jac, dotProd, t = _solution->getCurrentTime();
  double graduh[3]  = {0., 0., 0.};
  std::vector<double> gradu(3, 0.);

  if(error && _vectorSolution == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute H1 seminorm of error function"
                                " because exact solution is NULL");
    exit(-1);
  }

  ElementTransformation T;

  for(int iElm = 0; iElm < _nElm; ++iElm) {

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    this->initializeLocalSolutionOnSpace(0, iElm);

    for(int k = 0; k < _nQuad; ++k) {
      jac = _J[_nQuad * iElm + k];
      _cnc->computeElementTransformation(_geoCoord, k, jac, T);

      _spaces[0]->interpolateFieldAtQuadNode_physicalGradient(_localSol[0], k, T, graduh);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      if(error) _vectorSolution->eval(t, _pos, gradu);

      graduh[0] -= gradu[0];
      graduh[1] -= gradu[1];
      graduh[2] -= gradu[2];
      dotProd = graduh[0]*graduh[0] + graduh[1]*graduh[1] + graduh[2]*graduh[2];

      res += dotProd * jac * _w[k];
    }
  }
  return sqrt(res);
}

double feNorm::computeH1SemiNormErrorEstimator()
{
  double res = 0.0, jac, dotProd, t = _solution->getCurrentTime();
  double graduh[3] = {0., 0., 0.};
  double graduRec[3] = {0., 0., 0.};
  std::vector<double> gradu(3, 0.);

  if(_rec == nullptr) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot compute H1 error estimate "
                                " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  ElementTransformation T;

  for(int iElm = 0; iElm < _nElm; ++iElm) {

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    this->initializeLocalSolutionOnSpace(0, iElm);

    for(int k = 0; k < _nQuad; ++k) {
      jac = _J[_nQuad * iElm + k];
      _cnc->computeElementTransformation(_geoCoord, k, jac, T);

      _spaces[0]->interpolateFieldAtQuadNode_physicalGradient(_localSol[0], k, T, graduh);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      _vectorSolution->eval(t, _pos, gradu);

      graduRec[0] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 0, iElm, k);
      graduRec[1] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 1, iElm, k);
      graduRec[2] = 0.;

      gradu[2] = 0.;

      graduRec[0] -= gradu[0];
      graduRec[1] -= gradu[1];
      graduRec[2] -= gradu[2];
      dotProd = graduRec[0]*graduRec[0] + graduRec[1]*graduRec[1] + graduRec[2]*graduRec[2];

      res += dotProd * jac * _w[k];
    }
  }
  return sqrt(res);
}

double feNorm::computeH1Norm(bool error)
{
  double L2 = this->computeLpNorm(2, error);
  double SemiH1 = this->computeH1SemiNorm(2, error);
  return sqrt(L2*L2 + SemiH1*SemiH1);
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
  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh) reduction(+ : res) schedule(dynamic)
  // #endif
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

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, u) reduction(+ : res) schedule(dynamic)
  // #endif
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
  return 0.;
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
