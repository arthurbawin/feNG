#include "feNorm.h"

thread_local feFunctionArguments THREAD_LOCAL_ARGS;
thread_local std::vector<double> THREAD_LOCAL_POS(3, 0.);
thread_local std::vector<double> THREAD_LOCAL_GEOCOORD(18, 0.);

feStatus createNorm(feNorm                            *&norm,
                    normType                            type,
                    const std::vector<const feSpace *> &spaces,
                    feSolution                         *sol,
                    const feFunction                         *scalarSolution,
                    const feVectorFunction                   *vectorSolution)
{
  for (const feSpace *fS : spaces)
  {
    if (fS == nullptr)
    {
      return feErrorMsg(FE_STATUS_ERROR,
                        "Null pointer in vector of FE spaces, maybe you forgot "
                        "to initialize it.");
    }
  }

  // Check that all interpolation spaces are defined on the same connectivity
  int cncGeoTag = spaces[0]->getCncGeoTag();
  for (const feSpace *fS : spaces)
  {
    if (fS->getCncGeoTag() != cncGeoTag)
      return feErrorMsg(
        FE_STATUS_ERROR,
        "Norm defined on more than one Physical group (connectivity):\n"
        "\t Finite element space for field \"%s\" is defined on connectivity "
        "\"%s\"\n"
        "\t Finite element space for field \"%s\" is defined on connectivity "
        "\"%s\"\n"
        "All finite element spaces in \"spaces\" should be defined on the same "
        "connectivity.\n",
        spaces[0]->getFieldID().data(),
        spaces[0]->getCncGeoID().data(),
        fS->getFieldID().data(),
        fS->getCncGeoID().data());
  }

  norm = new feNorm(type, spaces, sol, scalarSolution, vectorSolution);

  return FE_STATUS_OK;
}

feNorm::feNorm(normType                            type,
               const std::vector<const feSpace *> &spaces,
               feSolution                         *sol,
               const feFunction                         *scalarSolution,
               const feVectorFunction                   *vectorSolution)
  : _type(type)
  , _spaces(spaces)
  , _solution(sol)
  , _cnc(spaces[0]->getCncGeo())
  , _J(_cnc->getJacobians())
  , _scalarSolution(scalarSolution)
  , _vectorSolution(vectorSolution)
{
  _nQuad    = spaces[0]->getNumQuadPoints();
  _w        = spaces[0]->getQuadratureWeights();
  _geoSpace = _cnc->getFeSpace();
  _nElm     = _spaces[0]->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * _cnc->getNumVerticesPerElem());

  _adr.resize(spaces.size());
  _localSol.resize(spaces.size());
  _localSolDot.resize(spaces.size());
  for (size_t k = 0; k < spaces.size(); ++k)
  {
    _adr[k].resize(spaces[k]->getNumFunctions());
    _localSol[k].resize(spaces[k]->getNumFunctions());
    _localSolDot[k].resize(spaces[k]->getNumFunctions());
  }

  // Count number of scalar spaces
  _nScalarSpaces = 0;
  for(const auto &s: _spaces)
    if(!s->isVectorValued())
      ++_nScalarSpaces;

  _args.userScalarFields.resize(_nScalarSpaces);

  // Provide only a geometric connectivity, to compute integrals of source
  // terms, etc. but not of the finite element solution.
  _cncOnly = false;
}

feNorm::feNorm(feCncGeo         *cnc,
               const feFunction       *scalarSolution,
               const feVectorFunction *vectorSolution)
  : _cnc(cnc)
  , _J(_cnc->getJacobians())
  , _scalarSolution(scalarSolution)
  , _vectorSolution(vectorSolution)
{
  // FIXME: do not exit from constructor
  if (!cnc->jacobiansWereComputed())
  {
    feErrorMsg(
      FE_STATUS_ERROR,
      "Jacobian determinants were not computed on Physical Entity \"%s\". "
      "Define a finite element space on this entity first to compute areas and "
      "integrals.",
      cnc->getID().data());
    exit(-1);
  }

  _nQuad    = cnc->getFeSpace()->getNumQuadPoints();
  _w        = cnc->getFeSpace()->getQuadratureWeights();
  _geoSpace = cnc->getFeSpace();
  _nElm     = cnc->getNumElements();
  _pos.resize(3);
  _geoCoord.resize(3 * cnc->getNumVerticesPerElem());

  // Provide only a geometric connectivity, to compute integrals of source
  // terms, etc. but not of the finite element solution.
  _cncOnly = true;
}

double feNorm::compute(normType type)
{
  if (_cncOnly && !(type == AREA || type == INTEGRAL_F))
  {
    feErrorMsg(FE_STATUS_ERROR,
               "This norm was created with a geometric"
               " connectivity only. Cannot compute this norm without"
               " providing solution and finite element space.");
    exit(-1);
  }

  double res = 0.;
  switch (type)
  {
    case NONE:
      feWarning("No compute() method available for feNorm with type "
                "NONE"
                "");
      return 0.;
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
    case LINF_UDOT_ERROR:
      res = this->computeLInfNormTimeDerivative(true);
      break;
    case VECTOR_LINF_ERROR:
      res = this->computeVectorLInfNorm(true);
      break;
    case VECTOR_LINF_UDOT_ERROR:
      res = this->computeVectorLInfNormTimeDerivative(true);
      break;
    case SEMI_H1:
      res = this->computeH1SemiNorm(false);
      break;
    case SEMI_H1_ERROR:
      res = this->computeH1SemiNorm(true);
      break;
    case SEMI_H1_ERROR_ESTIMATE:
      res = this->computeH1SemiNormErrorEstimator();
      break;
    case VECTOR_SEMI_H1_ERROR:
      res = this->computeVectorH1SemiNorm(true);
      break;
    case H1:
      res = this->computeH1Norm(false);
      break;
    case H1_ERROR:
      res = this->computeH1Norm(true);
      break;
    // case VECTOR_H1_ERROR: : Requires giving two vector-valued fields (solution and gradient tensor)
    //   res = this->computeVectorH1Norm(true);
      break;
    case AREA:
      res = this->computeArea();
      break;
    case INTEGRAL:
      return this->computeIntegral();
    case INTEGRAL_SOL_DEPENDENT_F:
      return this->computeIntegralSolutionDependentFunction();
    case VECTOR_INTEGRAL_SOL_DEPENDENT_F:
      feErrorMsg(FE_STATUS_ERROR, "This integral is vector-valued. "
        "Call norm->compute(result) with result a vector instead.");
      exit(-1);
    case INTEGRAL_DT:
      return this->computeIntegralDt();
    case INTEGRAL_F:
      return this->computeIntegralUserFunction();
    case DOT_PRODUCT:
      return this->computeIntegralDotProduct();
    case PRESSURE_LIFT_FORCE:
      return this->computePressureLift();
    case PRESSURE_DRAG_FORCE:
      return this->computePressureDrag();
    case VISCOUS_LIFT_FORCE:
      return this->computeViscousLift();
    case VISCOUS_DRAG_FORCE:
      feErrorMsg(FE_STATUS_ERROR,
                 "Viscous drag computation requires providing reconstructed u "
                 "and v velocities."
                 " Call with e.g. : norm->computeViscousDrag(viscosity, &recU, "
                 "&recV); ");
      return 0.;
  }
  return res;
}

void feNorm::compute(std::vector<double> &result)
{
  switch(_type) {
    case VECTOR_INTEGRAL_SOL_DEPENDENT_F:
      computeVectorIntegralSolutionDependentFunction(result);
      break;
    default:
      feErrorMsg(FE_STATUS_ERROR, "This norm or integral does not have "
        "a compute(std::vector<double> &result) implemented.");
      exit(-1);
  }
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

void feNorm::initializeLocalSolutionTimeDerivativeOnSpace(int iSpace, int iElm)
{
  _spaces[iSpace]->initializeAddressingVector(iElm, _adr[iSpace]);
  _solution->getSolDotAtDOF(_adr[iSpace], _localSolDot[iSpace]);
}

void feNorm::assignFieldsToCallbackArguments(const int iNode)
{
  // Assign scalar- and vector-valued space to callback arguments
  // Currently limited to one vector space
  bool vectorWasSet = false;
  int indexScalar = 0;
  for(size_t iSpace = 0; iSpace < _spaces.size(); ++iSpace)
  {
    const feSpace *s = _spaces[iSpace];

    if (s->isVectorValued())
    {
      if (vectorWasSet)
      {
        // Handles only one field for now
        feErrorMsg(
          FE_STATUS_ERROR,
          "Trying to add more than one vector-valued field to function "
          "argument. Can only handle a single field for now.");
        exit(-1);
      }
      vectorWasSet = true;
      s->interpolateVectorFieldAtQuadNode(_localSol[iSpace],
                                          iNode,
                                          _args.v,
                                          s->getNumComponents());
    }
    else
    {
      const double scalar_val = s->interpolateFieldAtQuadNode(_localSol[iSpace], iNode);
      if(indexScalar == 0)
        _args.u = scalar_val;
      _args.userScalarFields[indexScalar] = scalar_val;
      indexScalar++;
    }
  }
}

feStatus feNorm::probeScalarField(const std::vector<double> &pos, double &res)
{
  feMesh *mesh = _spaces[0]->getMeshPtr();
  int     elm  = -1;
  double  r[3];
  bool    isFound =
    static_cast<feMesh2DP1 *>(mesh)->locateVertex(pos.data(), elm, r);
  if (!isFound)
    return feErrorMsg(
      FE_STATUS_ERROR,
      "Could not probe because point (%f,%f,%f) could not be found in the mesh",
      pos[0],
      pos[1],
      pos[2]);

  initializeLocalSolutionOnSpace(0, elm);
  res = _spaces[0]->interpolateField(_localSol[0], r);

  return FE_STATUS_OK;
}

double feNorm::computeLpNorm(int p, bool error)
{
  double res = 0.0, uh, u;
  _args.t    = _solution->getCurrentTime();

  if (error && _scalarSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp norm of error function"
               " because exact solution is NULL");
    exit(-1);
  }

  FILE *myfile = nullptr;
  if (_plotErrorToFile)
  {
    myfile = fopen(_errorPlotFileName.data(), "w");
    fprintf(myfile, "View \"LpNorm_errorOnElements\"{\n");
  }

  std::vector<double> errorOnElementsToThePowerP(_nElm, 0.);

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    double eLocPowP = 0.;
    for (int k = 0; k < _nQuad; ++k)
    {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      u = error ? _scalarSolution->eval(_args) : 0.0;

      eLocPowP += pow(fabs(u - uh), p) * _J[_nQuad * iElm + k] * _w[k];
    }

    errorOnElementsToThePowerP[iElm] = eLocPowP;
    res += eLocPowP;

    double eLoc = pow(eLocPowP, 1. / (double)p);

    if (_plotErrorToFile)
    {
      // Plot error
      _cnc->writeElementToPOS(myfile, _geoCoord, eLoc);
    }
  }

  if (_plotErrorToFile)
  {
    fprintf(myfile, "};\n");
    fclose(myfile);
  }

  if (_ignoreHighestError)
  {
    // Ignore the top X percents of error contributions
    std::sort(errorOnElementsToThePowerP.begin(),
              errorOnElementsToThePowerP.end());
    size_t numElmToIgnore =
      ceil(errorOnElementsToThePowerP.size() * _fractionToIgnore);
    res = 0.;
    for (size_t i = 0; i < errorOnElementsToThePowerP.size() - numElmToIgnore;
         ++i)
      res += errorOnElementsToThePowerP[i];
  }

  // Negative quadrature weights can sum to a very small but negative integral
  if (res < 0 && fabs(res) < 1e-14)
  {
    res = fabs(res);
  }

  return pow(res, 1. / (double)p);
}

double feNorm::computeLpNormOnElement(int p, bool error, int iElm)
{
  double uh, u;
  _args.t = _solution->getCurrentTime();

  if (error && _scalarSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp norm of error function"
               " because exact solution is NULL");
    exit(-1);
  }

  this->initializeLocalSolutionOnSpace(0, iElm);
  _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

  double eLocPowP = 0.;
  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int k = 0; k < _nQuad; ++k)
  {
    uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
    _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
    u = error ? _scalarSolution->eval(_args) : 0.0;

    eLocPowP += pow(fabs(u - uh), p) * _J[_nQuad * iElm + k] * _w[k];
  }

  return eLocPowP;
}

double area(const std::vector<double> &triCoord)
{
  // radius of circle circumscribing a triangle
  double dist[3], k = 0.0;
  for (int i = 0; i < 3; i++)
  {
    double x0 = triCoord[3 * i + 0];
    double y0 = triCoord[3 * i + 1];
    double x1 = triCoord[3 * ((i + 1) % 3) + 0];
    double y1 = triCoord[3 * ((i + 1) % 3) + 1];
    dist[i]   = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    k += 0.5 * dist[i];
  }
  const double area =
    std::sqrt(k * (k - dist[0]) * (k - dist[1]) * (k - dist[2]));
  return area;
}

double feNorm::computeSquaredErrorOnElement(int iElm)
{
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  this->initializeLocalSolutionOnSpace(0, iElm);
  _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

  double eLocPowP = 0.;
#if defined(HAVE_OMP)
#  pragma omp parallel
#endif
  {
    double u, uh;

#if defined(HAVE_OMP)
#  pragma omp for reduction(+ : eLocPowP)
#endif
    for (int k = 0; k < _nQuad; ++k)
    {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord,
                                                  k,
                                                  THREAD_LOCAL_ARGS.pos);
      u = _scalarSolution->eval(THREAD_LOCAL_ARGS);
      eLocPowP += (u - uh) * (u - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }

  return eLocPowP;
}

double
feNorm::computeSquaredErrorFromEstimatorOnElement(int  iElm,
                                                  bool useAverageEvaluations)
{
  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  // double t = _solution->getCurrentTime();

  this->initializeLocalSolutionOnSpace(0, iElm);
  _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

  double eLocPowP = 0.;
  // #if defined(HAVE_OMP)
  // #pragma omp parallel
  // #endif
  {
    double uRec, uh;

    // #if defined(HAVE_OMP)
    // #pragma omp for reduction(+:eLocPowP)
    // #endif
    for (int k = 0; k < _nQuad; ++k)
    {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord,
                                                  k,
                                                  THREAD_LOCAL_POS);

      // uRec = _rec->evaluateRecoveryAtQuadNode(PPR::RECOVERY, 0, iElm, k);
      // #pragma omp critical // Localization in the mesh is not yet thread safe
      // :(
      uRec = _rec->evaluateRecovery(PPR::RECOVERY,
                                    0,
                                    THREAD_LOCAL_POS.data(),
                                    useAverageEvaluations);

      eLocPowP += (uRec - uh) * (uRec - uh) * _J[_nQuad * iElm + k] * _w[k];
    }
  }

  return eLocPowP;
}

double feNorm::computeLpErrorEstimator(int p)
{
  double res = 0.0, uh, uRec;

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
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
  double res = 0.0, uExact, uRec;
  _args.t    = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      uExact = _scalarSolution->eval(_args);
      uRec   = _rec->evaluateRecoveryAtQuadNode(PPR::RECOVERY, 0, iElm, k);

      res += pow(fabs(uRec - uExact), p) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

double feNorm::computeLInfErrorExactVsEstimator()
{
  double res = -1., uExact, uRec;
  _args.t    = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      uExact = _scalarSolution->eval(_args);
      uRec   = _rec->evaluateRecoveryAtQuadNode(PPR::RECOVERY, 0, iElm, k);

      res = fmax(res, fabs(uRec - uExact));
    }
  }
  return res;
}

double feNorm::computeLInfErrorExactVsEstimatorAtVertices()
{
  double normMax      = 0.0, uRec, uExact;
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  const std::vector<int> &verticesConnectivity =
    _cnc->getVerticesConnectivity();

  // Compute max norm at mesh vertices at some distance L from the boundary
  // (Guo, Zhang, Zhao paper)
  int nVerticesPerElem = _cnc->getNumVerticesPerElem();
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

    for (int iVert = 0; iVert < nVerticesPerElem; ++iVert)
    {
      THREAD_LOCAL_ARGS.pos[0] = THREAD_LOCAL_GEOCOORD[3 * iVert + 0];
      THREAD_LOCAL_ARGS.pos[1] = THREAD_LOCAL_GEOCOORD[3 * iVert + 1];
      THREAD_LOCAL_ARGS.pos[2] = THREAD_LOCAL_GEOCOORD[3 * iVert + 2];

      // Check if vertex is at distance L from the boundary
      uExact     = _scalarSolution->eval(THREAD_LOCAL_ARGS);
      int vertex = verticesConnectivity[nVerticesPerElem * iElm + iVert];
      uRec       = _rec->evaluateRecoveryAtVertex(PPR::RECOVERY, 0, vertex);
      normMax    = fmax(normMax, fabs(uExact - uRec));
    }
  }

  feInfo("Computed norm of error in %f s", toc());

  return normMax;
}

thread_local std::vector<double> GRAD_U_EXACT(2, 0.);

double feNorm::computeSemiH1ErrorExactVsEstimator(bool   excludeBoundary,
                                                  double boundaryWidth,
                                                  double boundingBox[4])
{
  double res          = 0.0;
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

#if defined(HAVE_OMP)
#  pragma omp parallel
#endif
  {
    double normSquared, gradRec[2];

#if defined(HAVE_OMP)
#  pragma omp for reduction(+ : res)
#endif
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

      for (int k = 0; k < _nQuad; ++k)
      {
        _geoSpace->interpolateVectorFieldAtQuadNode(THREAD_LOCAL_GEOCOORD,
                                                    k,
                                                    THREAD_LOCAL_ARGS.pos);

        // Check if vertex is at distance L from the boundary
        // Hard-coded for [0,1] x [0,1] box
        double x     = THREAD_LOCAL_ARGS.pos[0];
        double y     = THREAD_LOCAL_ARGS.pos[1];
        bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                      y > boundingBox[2] + boundaryWidth &&
                      x < boundingBox[1] - boundaryWidth &&
                      y < boundingBox[3] - boundaryWidth);

        if (!excludeBoundary)
          inBox = true;

        if (inBox)
        {
          _vectorSolution->eval(THREAD_LOCAL_ARGS, GRAD_U_EXACT);
          gradRec[0] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 0, iElm, k);
          gradRec[1] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 1, iElm, k);

          normSquared =
            (GRAD_U_EXACT[0] - gradRec[0]) * (GRAD_U_EXACT[0] - gradRec[0]) +
            (GRAD_U_EXACT[1] - gradRec[1]) * (GRAD_U_EXACT[1] - gradRec[1]);

          res += normSquared * _J[_nQuad * iElm + k] * _w[k];
        }
      }
    }
  }

  feInfo("Computed norm of gradient error in %f s", toc());

  return sqrt(res);
}

double feNorm::computeSemiH1ErrorExactVsEstimator_Linf(bool   excludeBoundary,
                                                       double boundaryWidth,
                                                       double boundingBox[4])
{
  double normMax      = 0.0, gradRec[2];
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  const std::vector<int> &verticesConnectivity =
    _cnc->getVerticesConnectivity();

  // Compute max norm at mesh vertices at some distance L from the boundary
  // (Guo, Zhang, Zhao paper)
  int nVerticesPerElem = _cnc->getNumVerticesPerElem();
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

    for (int iVert = 0; iVert < nVerticesPerElem; ++iVert)
    {
      THREAD_LOCAL_ARGS.pos[0] = THREAD_LOCAL_GEOCOORD[3 * iVert + 0];
      THREAD_LOCAL_ARGS.pos[1] = THREAD_LOCAL_GEOCOORD[3 * iVert + 1];
      THREAD_LOCAL_ARGS.pos[2] = THREAD_LOCAL_GEOCOORD[3 * iVert + 2];

      // Check if vertex is at distance L from the boundary
      // Hard-coded for [0,1] x [0,1] box
      double x     = THREAD_LOCAL_ARGS.pos[0];
      double y     = THREAD_LOCAL_ARGS.pos[1];
      bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                    y > boundingBox[2] + boundaryWidth &&
                    x < boundingBox[1] - boundaryWidth &&
                    y < boundingBox[3] - boundaryWidth);

      if (!excludeBoundary)
        inBox = true;

      if (inBox)
      {
        _vectorSolution->eval(THREAD_LOCAL_ARGS, GRAD_U_EXACT);

        int vertex = verticesConnectivity[nVerticesPerElem * iElm + iVert];

        gradRec[0] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
        gradRec[1] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

        for (int i = 0; i < 2; ++i)
        {
          normMax = fmax(normMax, fabs(GRAD_U_EXACT[i] - gradRec[i]));
        }
      }
    }
  }

  feInfo("Computed norm of gradient error in %f s", toc());

  return normMax;
}

thread_local std::vector<double> HESS_U_EXACT(4, 0.);

double feNorm::computeErrorHessianExactVsEstimator(bool   excludeBoundary,
                                                   double boundaryWidth,
                                                   double boundingBox[4])
{
  double res          = 0.;
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

#if defined(HAVE_OMP)
#  pragma omp parallel
#endif
  {
    double normSquared, hessRec[4];

#if defined(HAVE_OMP)
#  pragma omp for reduction(+ : res)
#endif
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

      for (int k = 0; k < _nQuad; ++k)
      {
        _geoSpace->interpolateVectorFieldAtQuadNode(THREAD_LOCAL_GEOCOORD,
                                                    k,
                                                    THREAD_LOCAL_ARGS.pos);

        // Check if vertex is at distance L from the boundary
        // Hard-coded for [0,1] x [0,1] box
        double x     = THREAD_LOCAL_ARGS.pos[0];
        double y     = THREAD_LOCAL_ARGS.pos[1];
        bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                      y > boundingBox[2] + boundaryWidth &&
                      x < boundingBox[1] - boundaryWidth &&
                      y < boundingBox[3] - boundaryWidth);

        if (!excludeBoundary)
          inBox = true;

        if (inBox)
        {
          _vectorSolution->eval(THREAD_LOCAL_ARGS, HESS_U_EXACT);

          // Averaged derivatives-wise
          hessRec[0] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 2, iElm, k);
          hessRec[1] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 3, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 4, iElm, k)) /
            2.;
          hessRec[2] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 5, iElm, k);

          normSquared =
            (HESS_U_EXACT[0] - hessRec[0]) * (HESS_U_EXACT[0] - hessRec[0]) +
            +2. * (HESS_U_EXACT[1] - hessRec[1]) *
              (HESS_U_EXACT[1] - hessRec[1]) +
            +(HESS_U_EXACT[3] - hessRec[2]) * (HESS_U_EXACT[3] - hessRec[2]);

          // Component-wise
          // hessRec[0] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 2,
          // iElm, k); hessRec[1] =
          // _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 3, iElm, k);
          // hessRec[2] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 4,
          // iElm, k); hessRec[3] =
          // _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 5, iElm, k);

          // normSquared = 0.;
          // for(int i = 0; i < 4; ++i)
          //   normSquared += (HESS_U_EXACT[i] - hessRec[i]) * (HESS_U_EXACT[i]
          //   - hessRec[i]);

          res += normSquared * _J[_nQuad * iElm + k] * _w[k];
        }
      }
    }
  }

  feInfo("Computed norm of hessian error in %f s", toc());

  return sqrt(res);
}

double feNorm::computeErrorHessianExactVsEstimator_Linf(bool   excludeBoundary,
                                                        double boundaryWidth,
                                                        double boundingBox[4])
{
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  double normMax = 0., hessRec[4];

  // Compute max norm at quadrature nodes
  // for(int iElm = 0; iElm < _nElm; ++iElm) {
  //   _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

  //   for(int k = 0; k < _nQuad; ++k) {
  //     _geoSpace->interpolateVectorFieldAtQuadNode(THREAD_LOCAL_GEOCOORD, k,
  //     THREAD_LOCAL_POS); _vectorSolution->eval(t, THREAD_LOCAL_POS,
  //     HESS_U_EXACT);

  //     // Component-wise
  //     hessRec[0] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 2, iElm,
  //     k); hessRec[1] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 3,
  //     iElm, k); hessRec[2] =
  //     _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 4, iElm, k);
  //     hessRec[3] = _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 5, iElm,
  //     k);

  //     for(int i = 0; i < 4; ++i) {
  //       normMax = fmax(normMax, fabs(HESS_U_EXACT[i] - hessRec[i]));
  //     }
  //   }
  // }

  const std::vector<int> &verticesConnectivity =
    _cnc->getVerticesConnectivity();

  // Compute max norm at mesh vertices at some distance L from the boundary
  // (Guo, Zhang, Zhao paper)
  int nVerticesPerElem = _cnc->getNumVerticesPerElem();
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

    for (int iVert = 0; iVert < nVerticesPerElem; ++iVert)
    {
      THREAD_LOCAL_ARGS.pos[0] = THREAD_LOCAL_GEOCOORD[3 * iVert + 0];
      THREAD_LOCAL_ARGS.pos[1] = THREAD_LOCAL_GEOCOORD[3 * iVert + 1];
      THREAD_LOCAL_ARGS.pos[2] = THREAD_LOCAL_GEOCOORD[3 * iVert + 2];

      // Check if vertex is at distance L from the boundary
      // Hard-coded for [0,1] x [0,1] box
      double x     = THREAD_LOCAL_ARGS.pos[0];
      double y     = THREAD_LOCAL_ARGS.pos[1];
      bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                    y > boundingBox[2] + boundaryWidth &&
                    x < boundingBox[1] - boundaryWidth &&
                    y < boundingBox[3] - boundaryWidth);

      if (!excludeBoundary)
        inBox = true;

      if (inBox)
      {
        _vectorSolution->eval(THREAD_LOCAL_ARGS, HESS_U_EXACT);

        int vertex = verticesConnectivity[nVerticesPerElem * iElm + iVert];

        hessRec[0] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
        hessRec[1] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex);
        hessRec[2] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex);
        hessRec[3] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

        for (int i = 0; i < 4; ++i)
        {
          normMax = fmax(normMax, fabs(HESS_U_EXACT[i] - hessRec[i]));
        }
      }
    }
  }

  feInfo("Computed norm of hessian error in %f s", toc());

  return normMax;
}

thread_local std::vector<double> D3U_EXACT(8, 0.);

double
feNorm::computeErrorThirdDerivativesExactVsEstimator(bool   excludeBoundary,
                                                     double boundaryWidth,
                                                     double boundingBox[4])
{
  double res          = 0.0;
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

#if defined(HAVE_OMP)
#  pragma omp parallel
#endif
  {
    double normSquared, d3Rec[4];

#if defined(HAVE_OMP)
#  pragma omp for reduction(+ : res)
#endif
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

      for (int k = 0; k < _nQuad; ++k)
      {
        _geoSpace->interpolateVectorFieldAtQuadNode(THREAD_LOCAL_GEOCOORD,
                                                    k,
                                                    THREAD_LOCAL_ARGS.pos);

        // Check if vertex is at distance L from the boundary
        // Hard-coded for [0,1] x [0,1] box
        double x     = THREAD_LOCAL_ARGS.pos[0];
        double y     = THREAD_LOCAL_ARGS.pos[1];
        bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                      y > boundingBox[2] + boundaryWidth &&
                      x < boundingBox[1] - boundaryWidth &&
                      y < boundingBox[3] - boundaryWidth);

        if (!excludeBoundary)
          inBox = true;

        if (inBox)
        {
          _vectorSolution->eval(THREAD_LOCAL_ARGS, D3U_EXACT);

          // Averaged derivatives-wise
          // Be careful for the order of mixed derivatives:
          // it's xxx, [xxy, xyx, xyy, yxx, yxy, yyx], yyy
          // so 7-8-10 and 9-11-12 (first 6 positions are 2 first and 4 second
          // derivatives)
          d3Rec[0] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 6, iElm, k);
          d3Rec[1] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 7, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 8, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 10, iElm, k)) /
            3.;
          d3Rec[2] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 9, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 11, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 12, iElm, k)) /
            3.;
          d3Rec[3] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 13, iElm, k);

          normSquared =
            (D3U_EXACT[0] - d3Rec[0]) * (D3U_EXACT[0] - d3Rec[0]) +
            3. * (D3U_EXACT[1] - d3Rec[1]) * (D3U_EXACT[1] - d3Rec[1]) +
            3. * (D3U_EXACT[3] - d3Rec[2]) * (D3U_EXACT[3] - d3Rec[2]) +
            (D3U_EXACT[7] - d3Rec[3]) * (D3U_EXACT[7] - d3Rec[3]);

          // // Component-wise
          // normSquared = 0.;
          // for(int ii = 0; ii < 8; ++ii) {
          //   d3Rec[ii] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE,
          //   6+ii, iElm, k); normSquared += (D3U_EXACT[ii] - d3Rec[ii]) *
          //   (D3U_EXACT[ii] - d3Rec[ii]);
          // }

          res += normSquared * _J[_nQuad * iElm + k] * _w[k];
        }
      }
    }
  }

  feInfo("Computed norm of d3u error in %f s", toc());

  return sqrt(res);
}

double
feNorm::computeErrorThirdDerivativesExactVsEstimator_Linf(bool excludeBoundary,
                                                          double boundaryWidth,
                                                          double boundingBox[4])
{
  double normMax      = 0.0, d3Rec[8];
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  const std::vector<int> &verticesConnectivity =
    _cnc->getVerticesConnectivity();

  // Compute max norm at mesh vertices at some distance L from the boundary
  // (Guo, Zhang, Zhao paper)
  int nVerticesPerElem = _cnc->getNumVerticesPerElem();
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

    for (int iVert = 0; iVert < nVerticesPerElem; ++iVert)
    {
      THREAD_LOCAL_ARGS.pos[0] = THREAD_LOCAL_GEOCOORD[3 * iVert + 0];
      THREAD_LOCAL_ARGS.pos[1] = THREAD_LOCAL_GEOCOORD[3 * iVert + 1];
      THREAD_LOCAL_ARGS.pos[2] = THREAD_LOCAL_GEOCOORD[3 * iVert + 2];

      // Check if vertex is at distance L from the boundary
      // Hard-coded for [0,1] x [0,1] box
      double x     = THREAD_LOCAL_ARGS.pos[0];
      double y     = THREAD_LOCAL_ARGS.pos[1];
      bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                    y > boundingBox[2] + boundaryWidth &&
                    x < boundingBox[1] - boundaryWidth &&
                    y < boundingBox[3] - boundaryWidth);

      if (!excludeBoundary)
        inBox = true;

      if (inBox)
      {
        _vectorSolution->eval(THREAD_LOCAL_ARGS, D3U_EXACT);

        int vertex = verticesConnectivity[nVerticesPerElem * iElm + iVert];

        d3Rec[0] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
        d3Rec[1] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex);
        d3Rec[2] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex);
        d3Rec[3] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex);
        d3Rec[4] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex);
        d3Rec[5] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex);
        d3Rec[6] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex);
        d3Rec[7] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);

        for (int i = 0; i < 8; ++i)
        {
          normMax = fmax(normMax, fabs(D3U_EXACT[i] - d3Rec[i]));
        }
      }
    }
  }

  feInfo("Computed norm of d3u error in %f s", toc());

  return normMax;
}

thread_local std::vector<double> D4U_EXACT(16, 0.);

double
feNorm::computeErrorFourthDerivativesExactVsEstimator(bool   excludeBoundary,
                                                      double boundaryWidth,
                                                      double boundingBox[4])
{
  double res          = 0.0;
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

#if defined(HAVE_OMP)
#  pragma omp parallel
#endif
  {
    double normSquared, d4Rec[5];

#if defined(HAVE_OMP)
#  pragma omp for reduction(+ : res)
#endif
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

      for (int k = 0; k < _nQuad; ++k)
      {
        _geoSpace->interpolateVectorFieldAtQuadNode(THREAD_LOCAL_GEOCOORD,
                                                    k,
                                                    THREAD_LOCAL_ARGS.pos);

        // Check if vertex is at distance L from the boundary
        // Hard-coded for [0,1] x [0,1] box
        double x     = THREAD_LOCAL_ARGS.pos[0];
        double y     = THREAD_LOCAL_ARGS.pos[1];
        bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                      y > boundingBox[2] + boundaryWidth &&
                      x < boundingBox[1] - boundaryWidth &&
                      y < boundingBox[3] - boundaryWidth);

        if (!excludeBoundary)
          inBox = true;

        if (inBox)
        {
          _vectorSolution->eval(THREAD_LOCAL_ARGS, D4U_EXACT);

          // Averaged derivatives-wise
          d4Rec[0] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 14, iElm, k);
          d4Rec[1] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 15, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 16, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 18, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 22, iElm, k)) /
            4.;
          d4Rec[2] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 17, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 19, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 20, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 23, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 24, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 26, iElm, k)) /
            6.;
          d4Rec[3] =
            (_rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 21, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 25, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 27, iElm, k) +
             _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 28, iElm, k)) /
            4.;
          d4Rec[4] =
            _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 29, iElm, k);

          normSquared =
            (D4U_EXACT[0] - d4Rec[0]) * (D4U_EXACT[0] - d4Rec[0]) +
            4. * (D4U_EXACT[1] - d4Rec[1]) * (D4U_EXACT[1] - d4Rec[1]) +
            6. * (D4U_EXACT[3] - d4Rec[2]) * (D4U_EXACT[3] - d4Rec[2]) +
            4. * (D4U_EXACT[7] - d4Rec[3]) * (D4U_EXACT[7] - d4Rec[3]) +
            (D4U_EXACT[15] - d4Rec[4]) * (D4U_EXACT[15] - d4Rec[4]);

          // // Component-wise
          // normSquared = 0.;
          // for(int ii = 0; ii < 8; ++ii) {
          //   d4Rec[ii] =  _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE,
          //   6+ii, iElm, k); normSquared += (D4U_EXACT[ii] - d4Rec[ii]) *
          //   (D4U_EXACT[ii] - d4Rec[ii]);
          // }

          res += normSquared * _J[_nQuad * iElm + k] * _w[k];
        }
      }
    }
  }

  feInfo("Computed norm of d4u error in %f s", toc());

  return sqrt(res);
}

double feNorm::computeErrorFourthDerivativesExactVsEstimator_Linf(
  bool   excludeBoundary,
  double boundaryWidth,
  double boundingBox[4])
{
  double normMax      = 0.0, d4Rec[16];
  THREAD_LOCAL_ARGS.t = _solution->getCurrentTime();

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  tic();

  const std::vector<int> &verticesConnectivity =
    _cnc->getVerticesConnectivity();

  FILE *file = fopen("testInside.pos", "w");
  fprintf(file, "View\"testInside\"{\n");

  FILE *myfile;
  myfile = fopen("errorFourthDerivative_Linf.pos", "w");
  fprintf(myfile, "View \"errorOnElements\"{\n");

  // Compute max norm at mesh vertices at some distance L from the boundary
  // (Guo, Zhang, Zhao paper)
  int nVerticesPerElem = _cnc->getNumVerticesPerElem();
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    double maxErrorOnElement = 0.;
    _spaces[0]->_mesh->getCoord(_cnc, iElm, THREAD_LOCAL_GEOCOORD);

    for (int iVert = 0; iVert < nVerticesPerElem; ++iVert)
    {
      THREAD_LOCAL_ARGS.pos[0] = THREAD_LOCAL_GEOCOORD[3 * iVert + 0];
      THREAD_LOCAL_ARGS.pos[1] = THREAD_LOCAL_GEOCOORD[3 * iVert + 1];
      THREAD_LOCAL_ARGS.pos[2] = THREAD_LOCAL_GEOCOORD[3 * iVert + 2];

      // Check if vertex is at distance L from the boundary
      // Hard-coded for [0,1] x [0,1] box
      double x     = THREAD_LOCAL_ARGS.pos[0];
      double y     = THREAD_LOCAL_ARGS.pos[1];
      bool   inBox = (x > boundingBox[0] + boundaryWidth &&
                    y > boundingBox[2] + boundaryWidth &&
                    x < boundingBox[1] - boundaryWidth &&
                    y < boundingBox[3] - boundaryWidth);

      if (!excludeBoundary)
        inBox = true;

      if (inBox)
      {
        fprintf(file, "SP(%g,%g,0.){%g};\n", x, y, 1.);
        _vectorSolution->eval(THREAD_LOCAL_ARGS, D4U_EXACT);

        int vertex = verticesConnectivity[nVerticesPerElem * iElm + iVert];

        d4Rec[0]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 14, vertex);
        d4Rec[1]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 15, vertex);
        d4Rec[2]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 16, vertex);
        d4Rec[3]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 17, vertex);
        d4Rec[4]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 18, vertex);
        d4Rec[5]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 19, vertex);
        d4Rec[6]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 20, vertex);
        d4Rec[7]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 21, vertex);
        d4Rec[8]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 22, vertex);
        d4Rec[9]  = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 23, vertex);
        d4Rec[10] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 24, vertex);
        d4Rec[11] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 25, vertex);
        d4Rec[12] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 26, vertex);
        d4Rec[13] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 27, vertex);
        d4Rec[14] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 28, vertex);
        d4Rec[15] = _rec->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 29, vertex);

        for (int i = 0; i < 16; ++i)
        {
          normMax = fmax(normMax, fabs(D4U_EXACT[i] - d4Rec[i]));
          maxErrorOnElement =
            fmax(maxErrorOnElement, fabs(D4U_EXACT[i] - d4Rec[i]));
        }

        // Plot error
        _cnc->writeElementToPOS(myfile,
                                THREAD_LOCAL_GEOCOORD,
                                maxErrorOnElement);
      }
    }
  }
  fprintf(file, "};\n");
  fclose(file);
  fprintf(myfile, "};\n");
  fclose(myfile);

  feInfo("Computed norm of d4u error in %f s", toc());

  return normMax;
}

double feNorm::computeLpErrorFromTransferredSolution(int         p,
                                                     feSolution *otherSol)
{
  double res = 0.0, uh1, uh2;
  _localSol.resize(2);
  _localSol[0].resize(_spaces[0]->getNumFunctions());
  _localSol[1].resize(_spaces[0]->getNumFunctions());

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    // Initialize current and other solutions
    this->initializeLocalSolutionOnSpace(0, iElm);
    otherSol->getSolAtDOF(_adr[0], _localSol[1]);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      uh1 = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      uh2 = _spaces[0]->interpolateFieldAtQuadNode(_localSol[1], k);

      res += pow(fabs(uh1 - uh2), p) * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

double feNorm::computeVectorLpNorm(int p, bool error)
{
  double res = 0.0, compNorm;
  _args.t    = _solution->getCurrentTime();
  std::vector<double> u(3, 0.);
  std::vector<double> uh(3, 0.);
  int                 nComponents = _spaces[0]->getNumComponents();

  if (error && _vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Lp norm of vector error function"
               " because exact solution is NULL");
    exit(-1);
  }

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      _spaces[0]->interpolateVectorFieldAtQuadNode(_localSol[0],
                                                   k,
                                                   uh,
                                                   nComponents);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      if (error)
        (*_vectorSolution)(_args, u);

      compNorm = 0.;
      for (int i = 0; i < nComponents; ++i)
      {
        compNorm += pow(fabs(u[i] - uh[i]), p);
      }

      res += compNorm * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return pow(res, 1. / (double)p);
}

double feNorm::computeVectorLInfNorm(bool error)
{
  double res = 0.0;
  _args.t    = _solution->getCurrentTime();
  std::vector<double> u(3, 0.);
  std::vector<double> uh(3, 0.);
  int                 nComponents = _spaces[0]->getNumComponents();

  if (error && _vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Linf norm of vector error function"
               " because exact solution is NULL");
    exit(-1);
  }

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      _spaces[0]->interpolateVectorFieldAtQuadNode(_localSol[0],
                                                   k,
                                                   uh,
                                                   nComponents);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      if (error)
        (*_vectorSolution)(_args, u);

      for (int i = 0; i < nComponents; ++i)
      {
        res = fmax(res, fabs(u[i] - uh[i]));
      }
    }
  }
  return res;
}

double feNorm::computeVectorLInfNormTimeDerivative(bool error)
{
  double res = 0.0;
  _args.t    = _solution->getCurrentTime();
  std::vector<double> dudt(3, 0.);
  std::vector<double> duhdt(3, 0.);
  int                 nComponents = _spaces[0]->getNumComponents();

  if (error && _vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute Linf norm of vector error function"
               " because exact solution is NULL");
    exit(-1);
  }

  double avg = 0.;
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionTimeDerivativeOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    //////////////////
    // const std::vector<double> &Lcoor = _spaces[0]->getLcoor();
    // for(size_t kk = 0; kk < Lcoor.size()/3; ++kk)
    // {
    //   // feInfo("%d", Lcoor.size());
    //   // feInfo("%d", _localSolDot[0].size());
    //   std::vector<double> xsi = {Lcoor[3*kk], Lcoor[3*kk+1], Lcoor[3*kk+2]};
    //   _geoSpace->interpolateVectorField(_geoCoord, xsi.data(), _args.pos);
    //   (*_vectorSolution)(_args, dudt);
    //   // _spaces[0]->interpolateVectorField(_localSolDot[0], k, duhdt,
    //   nComponents); feInfo("Exact dudt[0] at vertex %d = %f - dsoldt = %f",
    //   kk, dudt[kk % 2], _localSolDot[0][kk]);
    // }
    //////////////////

    for (int k = 0; k < _nQuad; ++k)
    {
      _spaces[0]->interpolateVectorFieldAtQuadNode(_localSolDot[0],
                                                   k,
                                                   duhdt,
                                                   nComponents);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      if (error)
        (*_vectorSolution)(_args, dudt);

      for (int i = 0; i < nComponents; ++i)
      {
        // feInfo("Comp %d = %f vs %f", i, dudt[i], duhdt[i]);
        res = fmax(res, fabs(dudt[i] - duhdt[i]));
        avg += fabs(dudt[i] - duhdt[i]) * _J[_nQuad * iElm + k] * _w[k];
      }
    }
  }
  // feInfo("Avg = %+-1.8e", avg);
  // feInfo("Err = %+-1.8e", res);
  return res;
}

// Deprecated by general Lp norm function
// double feNorm::computeL1Norm(bool error)
// {
//   double res = 0.0, uh, u, t = _solution->getCurrentTime();

//   // #if defined(HAVE_OMP)
//   // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res)
//   schedule(dynamic)
//   // #endif
//   for(int iElm = 0; iElm < _nElm; ++iElm) {
//     this->initializeLocalSolutionOnSpace(0, iElm);

//     _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

//     for(int k = 0; k < _nQuad; ++k) {
//       uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
//       _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
//       u = error ? _scalarSolution->eval(t, _pos) : 0.0;
//       res += fabs(u - uh) * _J[_nQuad * iElm + k] * _w[k];
//     }
//   }
//   return res;
// }

// Deprecated by general Lp norm function
// double feNorm::computeL2Norm(bool error)
// {
//   double res = 0.0, uh, u, t = _solution->getCurrentTime();

//   // #if defined(HAVE_OMP)
//   // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res)
//   schedule(dynamic)
//   // #endif
//   for(int iElm = 0; iElm < _nElm; ++iElm) {
//     this->initializeLocalSolutionOnSpace(0, iElm);

//     _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

//     for(int k = 0; k < _nQuad; ++k) {
//       uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
//       _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
//       u = error ? _scalarSolution->eval(t, _pos) : 0.0;
//       res += (u - uh) * (u - uh) * _J[_nQuad * iElm + k] * _w[k];
//     }
//   }
//   return sqrt(res);
// }

double feNorm::computeLInfNorm(bool error)
{
  double res = -DBL_MAX, uh, u;
  _args.t    = _solution->getCurrentTime();

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      u   = error ? _scalarSolution->eval(_args) : 0.0;
      res = fmax(res, fabs(u - uh));
    }
  }
  return res;
}

double feNorm::computeLInfNormTimeDerivative(bool error)
{
  double res = -DBL_MAX, duhdt, dudt;
  _args.t    = _solution->getCurrentTime();

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh, u) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionTimeDerivativeOnSpace(0, iElm);

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      duhdt = _spaces[0]->interpolateFieldAtQuadNode(_localSolDot[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      dudt = error ? _scalarSolution->eval(_args) : 0.0;
      res  = fmax(res, fabs(dudt - duhdt));
    }
  }
  return res;
}

double feNorm::computeH1SemiNorm(bool error)
{
  double              res       = 0.0, jac, dotProd;
  double              graduh[3] = {0., 0., 0.};
  std::vector<double> gradu(3, 0.);
  _args.t = _solution->getCurrentTime();

  if (error && _vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute H1 seminorm of error function"
               " because exact solution is NULL");
    exit(-1);
  }

  ElementTransformation T;

  FILE *myfile = nullptr;
  if (_plotErrorToFile)
  {
    myfile = fopen(_errorPlotFileName.data(), "w");
    fprintf(myfile, "View \"H1_seminorm_errorOnElements\"{\n");
  }

  std::vector<double> errorSquaredOnElements(_nElm, 0.);

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    if (_cnc->hasConstantTransformation())
    {
      _cnc->getElementTransformation(iElm, T);
    }

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    this->initializeLocalSolutionOnSpace(0, iElm);

    double eLocSquared = 0.;
    for (int k = 0; k < _nQuad; ++k)
    {
      jac = _J[_nQuad * iElm + k];
      _cnc->getElementTransformation(_geoCoord, k, jac, T);

      _spaces[0]->interpolateFieldAtQuadNode_physicalGradient(_localSol[0],
                                                              k,
                                                              T,
                                                              graduh);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      if (error)
        _vectorSolution->eval(_args, gradu);

      graduh[0] -= gradu[0];
      graduh[1] -= gradu[1];
      graduh[2] -= gradu[2];
      dotProd =
        graduh[0] * graduh[0] + graduh[1] * graduh[1] + graduh[2] * graduh[2];

      eLocSquared += dotProd * jac * _w[k];
    }

    errorSquaredOnElements[iElm] = eLocSquared;
    res += eLocSquared;

    double eLoc = sqrt(eLocSquared);
    if (_plotErrorToFile)
    {
      // Plot error
      _cnc->writeElementToPOS(myfile, _geoCoord, eLoc);
    }
  }

  if (_plotErrorToFile)
  {
    fprintf(myfile, "};\n");
    fclose(myfile);
  }

  if (_ignoreHighestError)
  {
    // Ignore the top X percents of error contributions
    std::sort(errorSquaredOnElements.begin(), errorSquaredOnElements.end());
    size_t numElmToIgnore =
      ceil(errorSquaredOnElements.size() * _fractionToIgnore);
    res = 0.;
    for (size_t i = 0; i < errorSquaredOnElements.size() - numElmToIgnore; ++i)
      res += errorSquaredOnElements[i];
  }

  return sqrt(res);
}

double feNorm::computeH1SemiNormErrorEstimator()
{
  double res         = 0.0, jac, dotProd; //, t = _solution->getCurrentTime();
  double graduh[3]   = {0., 0., 0.};
  double graduRec[3] = {0., 0., 0.};
  // std::vector<double> gradu(3, 0.);

  if (_rec == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute H1 error estimate "
               " because feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  ElementTransformation T;

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    if (_cnc->hasConstantTransformation())
    {
      _cnc->getElementTransformation(iElm, T);
    }

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    this->initializeLocalSolutionOnSpace(0, iElm);

    for (int k = 0; k < _nQuad; ++k)
    {
      jac = _J[_nQuad * iElm + k];
      _cnc->getElementTransformation(_geoCoord, k, jac, T);

      _spaces[0]->interpolateFieldAtQuadNode_physicalGradient(_localSol[0],
                                                              k,
                                                              T,
                                                              graduh);
      // _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      // _vectorSolution->eval(t, _pos, gradu);

      graduRec[0] =
        _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 0, iElm, k);
      graduRec[1] =
        _rec->evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, 1, iElm, k);
      graduRec[2] = 0.;

      // gradu[2] = 0.;

      graduRec[0] -= graduh[0];
      graduRec[1] -= graduh[1];
      graduRec[2] -= graduh[2];
      dotProd = graduRec[0] * graduRec[0] + graduRec[1] * graduRec[1] +
                graduRec[2] * graduRec[2];

      res += dotProd * jac * _w[k];
    }
  }
  return sqrt(res);
}

double feNorm::computeVectorH1SemiNorm(bool error)
{
  const int nComponents         = _spaces[0]->getNumComponents();
  std::vector<double> gradu(nComponents * nComponents, 0.);
  std::vector<double> graduh(nComponents * nComponents, 0.);
  double              res       = 0.0, jac, doubleContraction;
  _args.t = _solution->getCurrentTime();

  if (error && _vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute H1 seminorm because provided "
               "tensor-valued function for gradient is NULL");
    exit(-1);
  }

  ElementTransformation T;

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    if (_cnc->hasConstantTransformation())
    {
      _cnc->getElementTransformation(iElm, T);
    }

    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);
    this->initializeLocalSolutionOnSpace(0, iElm);

    for (int k = 0; k < _nQuad; ++k)
    {
      jac = _J[_nQuad * iElm + k];
      _cnc->getElementTransformation(_geoCoord, k, jac, T);
      _spaces[0]->interpolateVectorFieldAtQuadNode_physicalGradient(_localSol[0], nComponents, k, T, graduh.data());
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      if (error)
        _vectorSolution->eval(_args, gradu);

      // Convention : (grad(u))_ij := du_j/dx^i
      // 2D only for now (4 components in gradient tensor):
      graduh[0] -= gradu[0];
      graduh[1] -= gradu[1];
      graduh[2] -= gradu[2];
      graduh[3] -= gradu[3];
      doubleContraction =
          graduh[0] * graduh[0] + graduh[1] * graduh[1]
        + graduh[2] * graduh[2] + graduh[3] * graduh[3];

      res += doubleContraction * jac * _w[k];
    }
  }

  return sqrt(res);
}

double feNorm::computeH1Norm(bool error)
{
  double L2     = this->computeLpNorm(2, error);
  double SemiH1 = this->computeH1SemiNorm(error);
  return sqrt(L2 * L2 + SemiH1 * SemiH1);
}

// double feNorm::computeVectorH1Norm(bool error)
// {
//   double L2     = this->computeVectorLpNorm(2, error);
//   double SemiH1 = this->computeVectorH1SemiNorm(error);
//   return sqrt(L2 * L2 + SemiH1 * SemiH1);
// }

double feNorm::computeArea()
{
  double res = 0.0;
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    for (int k = 0; k < _nQuad; ++k)
    {
      res += _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeIntegral()
{
  double res = 0.0, uh;
  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, uh) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    for (int k = 0; k < _nQuad; ++k)
    {
      uh = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      res += uh * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeIntegralSolutionDependentFunction()
{
  if (_scalarSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute integral because vector-valued callback is NULL");
    exit(-1);
  }
  
  double res = 0., f;
  _args.t    = _solution->getCurrentTime();

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _geoSpace->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(size_t iSpace = 0; iSpace < _spaces.size(); ++iSpace)
      this->initializeLocalSolutionOnSpace(iSpace, iElm);

    for (int k = 0; k < _nQuad; ++k)
    {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      this->assignFieldsToCallbackArguments(k);

      // Evaluate the scalar-valued solution-dependent callback
      f = _scalarSolution->eval(_args);

      res += f * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

void feNorm::computeVectorIntegralSolutionDependentFunction(std::vector<double> &res)
{
  if (_vectorSolution == nullptr)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Cannot compute integral because vector-valued callback is NULL");
    exit(-1);
  }

  res[0] = res[1] = res[2] = 0.;
  double jacw;
  std::vector<double> F(3, 0.);
  _args.t = _solution->getCurrentTime();

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _geoSpace->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for(size_t iSpace = 0; iSpace < _spaces.size(); ++iSpace)
      this->initializeLocalSolutionOnSpace(iSpace, iElm);

    for (int k = 0; k < _nQuad; ++k)
    {
      jacw = _J[_nQuad * iElm + k] * _w[k];

      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      this->assignFieldsToCallbackArguments(k);

      // Evaluate the vector-valued solution-dependent callback
      _vectorSolution->eval(_args, F);

      res[0] += F[0] * jacw;
      res[1] += F[1] * jacw;
      res[2] += F[2] * jacw;
    }
  }
}

double feNorm::computeIntegralDt()
{
  double res = 0.0, uhDot;
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionTimeDerivativeOnSpace(0, iElm);
    for (int k = 0; k < _nQuad; ++k)
    {
      uhDot = _spaces[0]->interpolateFieldAtQuadNode(_localSolDot[0], k);
      res += uhDot * _J[_nQuad * iElm + k] * _w[k];
    }
  }
  return res;
}

double feNorm::computeIntegralUserFunction()
{
  double res = 0.0, u;
  _args.t    = 0.;

  // #if defined(HAVE_OMP)
  // #pragma omp parallel for private(_geoCoord, u) reduction(+ : res)
  // schedule(dynamic) #endif
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    _geoSpace->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _args.pos);
      u = (*_scalarSolution)(_args);
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
  //     _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(_geoCoord, k,
  //     dxdr);
  //     _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(_geoCoord, k,
  //     dxds); J = dxdr[0] * dxds[1] - dxdr[1] * dxds[0]; I +=
  //     _scalarSolution->eval(t, x) * solInt * J * _w[k];
  //   }
  // }
  // return I;
  return 0.;
}

/* Estimates the L2 norm of the error taking an external solution as the
 * reference solution. */
feStatus feNorm::computeErrorNormFromExternalSolution(
  feMesh                      *referenceMesh,
  const feSolution            *referenceSolution,
  const std::vector<feSpace *> referenceSpaces,
  feMesh                      *currentMesh,
  const feSolution            *currentSolution,
  std::vector<double>         &res)
{
  // Check that solutions time agree
  const double maxTimeDifference = 1e-10;
  if (fabs(currentSolution->getCurrentTime() -
           referenceSolution->getCurrentTime()) > maxTimeDifference)
  {
    return feErrorMsg(FE_STATUS_ERROR,
                      "Solutions time differ by more than %1.3e :\n"
                      "Time of other (reference) solution = %+-1.4e\n"
                      "Time of current           solution = %+-1.4e\n",
                      maxTimeDifference,
                      referenceSolution->getCurrentTime(),
                      currentSolution->getCurrentTime());
    // feWarning("Maybe the solution file does not match the current
    // solution.\n");
  }

  // Check that all spaces are found in referenceSpaces
  std::map<const feSpace *, const feSpace *> space2space;
  for (const auto &s : _spaces)
  {
    if (s->getNumComponents() > 1)
    {
      return feErrorMsg(FE_STATUS_ERROR, "Only computing scalar-valued error");
    }
    bool spaceWasFound = false;
    for (const auto &refs : referenceSpaces)
    {
      if (s->representsSameFieldAs(*refs))
      {
        spaceWasFound  = true;
        space2space[s] = refs;
        break;
      }
    }
    if (!spaceWasFound)
      return feErrorMsg(FE_STATUS_ERROR,
                        "Could not compute error from external solution!\n"
                        " Current finite element space for field \"%s\""
                        " on connectivity \"%s\" could not be found in the "
                        "given reference vector of spaces.",
                        s->getFieldID().data(),
                        s->getCncGeoID().data());
  }

  res.resize(_spaces.size(), 0.);
  std::vector<double> pos(3, 0.0);

  // Compute L2 error
  for (size_t iSpace = 0; iSpace < _spaces.size(); ++iSpace)
  {
    const feSpace *s        = _spaces[iSpace];
    const feSpace *refSpace = space2space.at(s);

    int    elm;
    double xsi[3], L2error = 0.0, uh, uhRef;

    std::vector<feInt>  refAddress(refSpace->getNumFunctions());
    std::vector<double> refLocalSol(refSpace->getNumFunctions());

    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      currentMesh->getCoord(_cnc, iElm, _geoCoord);

      // Initialize current solution on the element
      this->initializeLocalSolutionOnSpace(iSpace, iElm);

      for (int k = 0; k < _nQuad; ++k)
      {
        _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, pos);

        // Current solution
        uh = _spaces[iSpace]->interpolateFieldAtQuadNode(_localSol[iSpace], k);

        // Localize quadrature node in reference mesh
        // and evaluate reference solution
        bool isFound = static_cast<feMesh2DP1 *>(referenceMesh)
                         ->locateVertex(pos.data(), elm, xsi);

        if (!isFound)
        {
          feWarning("Point (%f, %f, %f) was not found in the reference mesh.\n",
                    pos[0],
                    pos[1],
                    pos[2]);
        }
        else
        {
          refSpace->initializeAddressingVector(elm, refAddress);
          referenceSolution->getSolAtDOF(refAddress, refLocalSol);

          // Reference solution at quadrature node
          uhRef = refSpace->interpolateField(refLocalSol, xsi);

          // Error
          L2error +=
            (uh - uhRef) * (uh - uhRef) * _J[_nQuad * iElm + k] * _w[k];
        }
      }
    }

    res[iSpace] = sqrt(L2error);
  }

  return FE_STATUS_OK;
}

// Compute lift force by evaluating integral of provided field (assumed to be
// pressure) projected along vertical unit vector.
double feNorm::computePressureLift()
{
  double res = 0.0, ph;

  if (_cnc->getDim() != 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute lift force on 1D connectivities for now.");
  }

  FILE *myf = fopen("NORMALES.pos", "w");
  fprintf(myf, "View \"normales\"{\n");

  std::vector<double> normalVectors; // Size = 3 * nQuad * nElm
  _cnc->computeNormalVectors(normalVectors);

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      ph = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);

      double ny = normalVectors[3 * _nQuad * iElm + 3 * k + 1];

      res += ph * ny * _J[_nQuad * iElm + k] * _w[k];
      fprintf(myf,
              "VP(%.16g,%.16g,%.16g){%.16g,%.16g,%.16g};\n",
              _pos[0],
              _pos[1],
              0.,
              normalVectors[3 * _nQuad * iElm + 3 * k + 0],
              normalVectors[3 * _nQuad * iElm + 3 * k + 1],
              normalVectors[3 * _nQuad * iElm + 3 * k + 2]);
    }
  }

  fprintf(myf, "};\n");
  fclose(myf);

  return res;
}

double feNorm::computeViscousLift()
{
  // To implement
  feErrorMsg(FE_STATUS_ERROR, "Viscous lift not implemented.");
  return 0.;
}

// Compute drag force by evaluating integral of provided field (assumed to be
// pressure) projected along horizontal unit vector.
double feNorm::computePressureDrag()
{
  double res = 0.0, ph;

  if (_cnc->getDim() != 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute drag force on 1D connectivities for now.");
  }

  std::vector<double> normalVectors; // Size = 3 * nQuad * nElm
  _cnc->computeNormalVectors(normalVectors);

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      ph = _spaces[0]->interpolateFieldAtQuadNode(_localSol[0], k);
      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);
      double nx = normalVectors[3 * _nQuad * iElm + 3 * k + 0];
      res += ph * nx * _J[_nQuad * iElm + k] * _w[k];
    }
  }

  return res;
}

double feNorm::computeViscousDrag(double         viscosity,
                                  feNewRecovery *recU,
                                  feNewRecovery *recV)
{
  double res = 0.0; //, uh, vh;

  if (_cnc->getDim() != 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute drag force on 1D connectivities for now.");
    exit(-1);
  }
  if (_spaces[0]->getNumComponents() == 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute drag force for a vector-valued FE space "
               "(velocity field).");
    exit(-1);
  }
  if (recU == nullptr || recV == nullptr)
  {
    feErrorMsg(
      FE_STATUS_ERROR,
      "Cannot computeViscousDrag "
      " because at least one feRecovery (solution reconstruction) is NULL");
    exit(-1);
  }

  std::vector<double> normalVectors; // Size = 3 * nQuad * nElm
  _cnc->computeNormalVectors(normalVectors);

  // double grad_u[4];
  // ElementTransformation transformation;

  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      // double jac = _J[_nQuad * iElm + k];
      // _cnc->computeElementTransformation(_geoCoord, k, jac, transformation);
      // _spaces[0]->interpolateVectorFieldAtQuadNode_physicalGradient(_localSol[0],
      // 2, k, transformation, grad_u); double dudx = grad_u[0]; double dudy =
      // grad_u[1]; double dvdx = grad_u[2]; double dvdy = grad_u[3];

      _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, k, _pos);

      double dudx = recU->evaluateRecovery(PPR::DERIVATIVE, 0, _pos.data());
      double dudy = recU->evaluateRecovery(PPR::DERIVATIVE, 1, _pos.data());
      double dvdx = recV->evaluateRecovery(PPR::DERIVATIVE, 0, _pos.data());
      double dvdy = recV->evaluateRecovery(PPR::DERIVATIVE, 1, _pos.data());

      double nx = normalVectors[3 * _nQuad * iElm + 3 * k + 0];
      double ny = normalVectors[3 * _nQuad * iElm + 3 * k + 1];

      double tauDotN[2] = {viscosity * (nx * (2. * dudx) + ny * (dudy + dvdx)),
                           viscosity * (nx * (dudy + dvdx) + ny * (2. * dvdy))};

      res += tauDotN[0] * nx * _J[_nQuad * iElm + k] * _w[k];
    }
  }

  return res;
}

double feNorm::computeForcesFromLagrangeMultiplier(int iComponent)
{
  double res = 0.0;

  if (_cnc->getDim() != 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute forces on 1D connectivities for now.");
    exit(-1);
  }
  if (_spaces[0]->getNumComponents() == 1)
  {
    feErrorMsg(FE_STATUS_ERROR,
               "Can only compute forces from a vector-valued FE space "
               "(lagrange multiplier).");
    exit(-1);
  }

  double lambda;
  for (int iElm = 0; iElm < _nElm; ++iElm)
  {
    this->initializeLocalSolutionOnSpace(0, iElm);
    _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

    for (int k = 0; k < _nQuad; ++k)
    {
      lambda = _spaces[0]->interpolateVectorFieldComponentAtQuadNode_fullField(
        _localSol[0], k, iComponent);
      res += lambda * _J[_nQuad * iElm + k] * _w[k];
    }
  }

  return res;
}

feStatus feNorm::reconstructScalarIsoline(const double val,
                                          std::vector<double> &isoline,
                                          double *userPtr)
{
  if (_spaces[0]->isVectorValued())
  {
    return feErrorMsg(FE_STATUS_ERROR, "Isoline reconstruction is only "
      "available for scalar fields");
  }

  isoline.clear();
  isoline.reserve(3 * 2 * _nElm); // Reserve for 2 intersection points per element (3 coord each)

  if(_cnc->getInterpolant() == geometricInterpolant::LINEP1)
  {
    std::vector<std::vector<double>> intersections;
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      this->initializeLocalSolutionOnSpace(0, iElm);
      _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

      // Check for sign changes in the nodal values
      bool signChanged = false;
      for(size_t i = 0; i < _localSol[0].size(); ++i) {
        for(size_t j = 0; j < _localSol[0].size(); ++j) {
          if( (_localSol[0][i] - val) * (_localSol[0][j] - val) < 0) {
            signChanged = true;
            break;
          }
        }
      }

      if(!signChanged) continue;

      intersections.reserve(1);

      // fprintf(elements, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n",
      //   _geoCoord[0], _geoCoord[1], _geoCoord[2],
      //   _geoCoord[3], _geoCoord[4], _geoCoord[5], 1., 1.);

      int numIntersections = 0;

      double xsi[3] = {0., 0., 0.};
      if(_spaces[0]->getRootOnEdge(0, _localSol[0], val, xsi, 1e-6))
      {
        // feInfo("Elm %d - Found root on edge 0 : (%+-1.4e, %+-1.4e, %+-1.4e)", iElm, xsi[0], xsi[1], xsi[2]);
        std::vector<double> pos(3);
        _geoSpace->interpolateVectorField(_geoCoord, xsi, pos);
        // fprintf(isolinePos, "SP(%g,%g,%g){%g};\n", pos[0], pos[1], pos[2], 2.);

        numIntersections++;
        intersections.push_back(pos);
        isoline.push_back(pos[0]);
        isoline.push_back(pos[1]);
        isoline.push_back(pos[2]);

        // Draw the intersection point
        // fprintf(isolinePos, "SP(%g,%g,%g){%g};\n",
        //   intersections[0][0], intersections[0][1], intersections[0][2], 1.);
      }

      intersections.clear();
    }
  }
  else if(_cnc->getInterpolant() == geometricInterpolant::TRIP1)
  {
    FILE *isolineTxt  = fopen("reconstructedIsoline.txt", "w");
    FILE *isolinePos  = fopen("reconstructedIsoline.pos", "w");
    FILE *elements = fopen("crossedElements.pos", "w");
    fprintf(isolinePos, "View \"reconstructedIsoline_val%f\"{\n", val);
    fprintf(elements, "View \"crossedElements\"{\n");

    std::vector<std::vector<double>> intersections;
    for (int iElm = 0; iElm < _nElm; ++iElm)
    {
      this->initializeLocalSolutionOnSpace(0, iElm);
      _spaces[0]->_mesh->getCoord(_cnc, iElm, _geoCoord);

      // Check for sign changes in the nodal values
      bool signChanged = false;
      for(size_t i = 0; i < _localSol[0].size(); ++i) {
        for(size_t j = 0; j < _localSol[0].size(); ++j) {
          if( (_localSol[0][i] - val) * (_localSol[0][j] - val) < 0) {
            signChanged = true;
            break;
          }
        }
      }

      if(!signChanged) continue;

      intersections.reserve(2);

      fprintf(elements, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n", _geoCoord[0],
                _geoCoord[1], _geoCoord[2], _geoCoord[3], _geoCoord[4], _geoCoord[5],
                _geoCoord[6], _geoCoord[7], _geoCoord[8], 1., 1., 1.);

      int numIntersections = 0;

      double xsi[3] = {0., 0., 0.};
      bool s0 = _spaces[0]->getRootOnEdge(0, _localSol[0], val, xsi, 1e-6);
      if(s0) {
        // feInfo("Elm %d - Found root on edge 0 : (%+-1.4e, %+-1.4e, %+-1.4e)", iElm, xsi[0], xsi[1], xsi[2]);
        std::vector<double> pos(3);
        _geoSpace->interpolateVectorField(_geoCoord, xsi, pos);
        fprintf(isolinePos, "SP(%g,%g,%g){%g};\n", pos[0], pos[1], pos[2], 2.);

        numIntersections++;
        intersections.push_back(pos);
        isoline.push_back(pos[0]);
        isoline.push_back(pos[1]);
        isoline.push_back(pos[2]);
      }

      bool s1 = _spaces[0]->getRootOnEdge(1, _localSol[0], val, xsi, 1e-6);
      if(s1) {
        // feInfo("Elm %d - Found root on edge 1 : (%+-1.4e, %+-1.4e, %+-1.4e)", iElm, xsi[0], xsi[1], xsi[2]);
        std::vector<double> pos(3);
        _geoSpace->interpolateVectorField(_geoCoord, xsi, pos);
        fprintf(isolinePos, "SP(%g,%g,%g){%g};\n", pos[0], pos[1], pos[2], 2.);

        numIntersections++;
        intersections.push_back(pos);
        isoline.push_back(pos[0]);
        isoline.push_back(pos[1]);
        isoline.push_back(pos[2]);
      }
      
      bool s2 = _spaces[0]->getRootOnEdge(2, _localSol[0], val, xsi, 1e-6);
      if(s2) {
        // feInfo("Elm %d - Found root on edge 2 : (%+-1.4e, %+-1.4e, %+-1.4e)", iElm, xsi[0], xsi[1], xsi[2]);
        std::vector<double> pos(3);
        _geoSpace->interpolateVectorField(_geoCoord, xsi, pos);
        fprintf(isolinePos, "SP(%g,%g,%g){%g};\n", pos[0], pos[1], pos[2], 2.);

        numIntersections++;
        intersections.push_back(pos);
        isoline.push_back(pos[0]);
        isoline.push_back(pos[1]);
        isoline.push_back(pos[2]);
      }

      if(numIntersections == 2) {
        // Isoline crosses the element once
        // Draw the intersection
        fprintf(isolinePos, "SL(%g,%g,%g,%g,%g,%g){%g,%g};\n",
          intersections[0][0], intersections[0][1], intersections[0][2],
          intersections[1][0], intersections[1][1], intersections[1][2], 1., 1.);

        // Print to text file in 2D
        fprintf(isolineTxt, "%+-1.16e \t %+-1.16e\n", intersections[0][0], intersections[0][1]);
        fprintf(isolineTxt, "%+-1.16e \t %+-1.16e\n", intersections[1][0], intersections[1][1]);

        // Get contact angle for particular case of Jurin's geometry
        // if(fabs(intersections[0][0] - 1.) < 1e-8 || fabs(intersections[1][0] - 1.) < 1e-8)
        if(fabs(intersections[0][1]) < 1e-8 || fabs(intersections[1][1]) < 1e-8)
        {
          double v1[2] = {intersections[1][0] - intersections[0][0],
                          intersections[1][1] - intersections[0][1]};
          const double normV1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1]);
          v1[0] /= normV1;
          v1[1] /= normV1;
          const double e_bdr[2] = {1.,  0.};
          // const double e_bdr[2] = {0., -1.};

          fprintf(isolinePos, "VP(%g,%g,%g){%g,%g,%g};\n", intersections[0][0], intersections[0][1], 0., v1[0], v1[1], 0.);
          fprintf(isolinePos, "VP(%g,%g,%g){%g,%g,%g};\n", intersections[0][0], intersections[0][1], 0., e_bdr[0], e_bdr[1], 0.);

          const double dotprod = v1[0]*e_bdr[0] + v1[1]*e_bdr[1];
          const double angle = acos(dotprod);
          if(userPtr != nullptr) {
            userPtr[0] = angle / M_PI * 180.;
          }
          feInfo("Angle is %f (deg)", angle / M_PI * 180.);
          feInfo("Angle is %f (deg)", angle / M_PI * 180.);
          feInfo("Angle is %f (deg)", angle / M_PI * 180.);
          feInfo("Angle is %f (deg)", angle / M_PI * 180.);
        }

      } else {
        feWarning("Isoline reconstruction: "
          "There are %d intersection(s) with this element", numIntersections);
        fprintf(elements, "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n", _geoCoord[0],
                _geoCoord[1], _geoCoord[2], _geoCoord[3], _geoCoord[4], _geoCoord[5],
                _geoCoord[6], _geoCoord[7], _geoCoord[8], 2., 2., 2.);
      }

      intersections.clear();
    }

    fprintf(isolinePos, "};");
    fprintf(elements, "};");
    fclose(isolineTxt);
    fclose(isolinePos);
    fclose(elements);
  }
  else
  {
    return feErrorMsg(FE_STATUS_ERROR, "Isoline extraction for cnc of type %s "
      "is not implemented", toString(_cnc->getInterpolant()).data());
  }

  return FE_STATUS_OK;
}