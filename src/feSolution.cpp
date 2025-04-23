#include "feSolution.h"

#include <iostream>
#include <fstream>

feSolution::feSolution(int numDOF, const std::vector<feSpace *> &spaces,
                       const std::vector<feSpace *> &essentialSpaces)
  : _nDOF(numDOF), _c0(0.), _tn(0.), _spaces(spaces), _essentialSpaces(essentialSpaces)
{
  _sol.resize(numDOF);
  _dsoldt.resize(numDOF);
  for(int i = 0; i < numDOF; ++i) {
    _sol[i] = 0.0;
    _dsoldt[i] = 0.0;
  }
}

feSolution::feSolution(const feSolutionContainer &container,
                       const int solutionIndex,
                       const std::vector<feSpace*> &spaces,
                       const std::vector<feSpace*> &essentialSpaces)
  : _nDOF(container.getNbDOFs()),
  _sol(container.getSolution(solutionIndex)),
  _dsoldt(container.getSolutionDot(solutionIndex)),
  _c0(0.), _tn(0.),
  _spaces(spaces), _essentialSpaces(essentialSpaces)
{

}

/* Constructs an feSolution from a file created by feSolution::printSol. */
feSolution::feSolution(std::string solutionFile)
  : _c0(0.), _spaces(std::vector<feSpace *>()), _essentialSpaces(std::vector<feSpace *>())
{
  feInfo("Reading solution file : %s\n", solutionFile.c_str());
  std::filebuf fb;
  if(fb.open(solutionFile, std::ios::in)) {
    std::istream input(&fb);
    std::string buffer;
    double solutionTime;
    input >> solutionTime;
    getline(input, buffer);
    input >> _nDOF;

    _tn = solutionTime;

    _sol.resize(_nDOF);
    _dsoldt.resize(_nDOF);

    double val;
    for(int i = 0; i < _nDOF; ++i) {
      getline(input, buffer);
      input >> val;
      _sol[i] = val;
    }
    for(int i = 0; i < _nDOF; ++i) {
      getline(input, buffer);
      input >> val;
      _dsoldt[i] = val;
    }

    fb.close();
  } else {
    feErrorMsg(FE_STATUS_ERROR, "Solution file could not be opened.\n");
    exit(-1);
  }
}

void feSolution::getSolAtDOF(const std::vector<feInt> &addressing, std::vector<double> &sol) const
{
  for(size_t i = 0; i < addressing.size(); ++i) {
    sol[i] = _sol[addressing[i]];
  }
}

void feSolution::getSolDotAtDOF(const std::vector<feInt> &addressing, std::vector<double> &solDot) const
{
  for(size_t i = 0; i < addressing.size(); ++i) {
    solDot[i] = _dsoldt[addressing[i]];
  }
}

void feSolution::initializeTemporalSolution(double t0, double t1, int nTimeSteps)
{
  _t0 = t0;
  _t1 = t1;
  _nTimeSteps = nTimeSteps;
  _dt = (t1 - t0) / (double)nTimeSteps;
  _tn = t0;
}

void feSolution::initialize(feMesh *mesh)
{
  this->initializeUnknowns(mesh);
  this->initializeEssentialBC(mesh);
}

void feSolution::initializeUnknowns(feMesh *mesh)
{
  for(feSpace *fS : _spaces) {
    if(fS->getDOFInitialization() == dofInitialization::PREVIOUS_SOL) {
      continue;
    }

    #if defined(HAVE_OMP)
    #pragma omp parallel
    #endif
    {
      // std::vector<double> x(3);
      feFunctionArguments args(_tn);
      std::vector<double> vecVal(3);
      double val;
      int nElm = fS->getNumElements();
      std::vector<feInt> adr(fS->getNumFunctions());
      std::vector<double> coor = fS->getLcoor();
      int numVerticesPerElem = mesh->getCncGeoByName(fS->getCncGeoID())->getNumVerticesPerElem();
      std::vector<double> localCoord(3 * numVerticesPerElem);

      feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());

      if(fS->getDOFInitialization() == dofInitialization::NODEWISE) {
        // Node based: the initial condition is imposed at the vertices DOF

        #if defined(HAVE_OMP)
        #pragma omp for
        #endif
        for(int iElm = 0; iElm < nElm; ++iElm) {
          fS->initializeAddressingVector(iElm, adr);
          mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

          for(int j = 0; j < fS->getNumFunctions(); ++j)
          {
            double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
            geoSpace->interpolateVectorField(localCoord, r, args.pos);

            if(fS->getNumComponents() > 1) {
              // Vector valued finite element space
              fS->evalFun(args, vecVal);
              int indexComponent = j % fS->getNumComponents();
              val = vecVal[indexComponent];
            } else {
              // Scalar valued finite element space
              val = fS->evalFun(args);
            }

            _sol[adr[j]] = val;

            /////////////////////////////////////////
            // // Solution time derivative
            // const double x = args.pos[0];
            // const double y = args.pos[1];
            // const double t = args.t;
            // const double a = 0.;
            // const double b = 0.;
            // vecVal[0] = -exp(-t) * (a * pow(x, 3) + x * y*y);
            // vecVal[1] = -exp(-t) * b * x*y;
            // vecVal[2] = 0.;

            // // const double t = args.t;
            // // vecVal[0] = -exp(-t);
            // // vecVal[1] = -exp(-t);
            // // vecVal[2] = 0.;

            // int indexComponent = j % fS->getNumComponents();
            // val = vecVal[indexComponent];
            // _dsoldt[adr[j]] = val;

            // Scalar valued time derivative
            // const double x = args.pos[0];
            // const double t = args.t;
            // const double udot = -2./((1.+t)*(1.+t)) * x;
            // const double udot = -2./((1.+t)*(1.+t));
            // const double udot = (double) rand() / (double)RAND_MAX;
            // const double udot = nan("");
            // _dsoldt[adr[j]] = udot;
            /////////////////////////////////////////

          }
        }
      }

      else if(fS->getDOFInitialization() == dofInitialization::LEAST_SQUARES) {
        // Assuming 1D Legendre polynomials: the initial condition is satisfied on each element in the
        // least-square sense:
        //
        //  <phi_i, phi_j> U_j^e = <phi_i, funSol> ==> A_ij^e U_j^e = B_j^e
        //
        // !!! We assume that the jacobian cancels out from both side, which is true only for P1 line
        // geometry !!!

        // The analytic diagonal mass matrix A_ii = 2/(2*i+1)
        std::vector<double> A(fS->getNumFunctions(), 0.);
        std::vector<double> B(fS->getNumFunctions(), 0.);
        for(int i = 0; i < fS->getNumFunctions(); ++i) {
          A[i] = 2. / (2. * i + 1.);
        }

        int nQuad = geoSpace->getNumQuadPoints();
        std::vector<double> wQuad = geoSpace->getQuadratureWeights();
        std::vector<double> xQuad = geoSpace->getRQuadraturePoints();

        double uQuad = 0.;
        std::vector<double> phi(fS->getNumFunctions(), 0.);
        for(int iElm = 0; iElm < nElm; ++iElm) {
          fS->initializeAddressingVector(iElm, adr);
          mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

          // Compute the RHS
          for(int i = 0; i < fS->getNumFunctions(); ++i) {
            B[i] = 0.;
            for(int k = 0; k < nQuad; ++k) {
              // Evaluate analytic solution at quad points
              double r[3] = {xQuad[k], 0., 0.};
              geoSpace->interpolateVectorField(localCoord, r, args.pos);
              uQuad = fS->evalFun(args);
              // Evaluate Legendre polynomials at quad points
              fS->L(r, phi.data());
              B[i] += wQuad[k] * phi[i] * uQuad;
            }
          }

          // Solve Aij Uj = Bj on the element
          for(int i = 0; i < fS->getNumFunctions(); ++i) {
            _sol[adr[i]] = B[i] / A[i];
          }
        }
      }
    } // pragma omp parallel
  }
}

void feSolution::initializeEssentialBC(feMesh *mesh, feSolutionContainer *solContainer)
{
  // std::vector<double> x(3);
  feFunctionArguments args(_tn);
  std::vector<double> vecVal(3);
  double val;

  // Essential spaces + non essential vector spaces with one or more
  // essential component. The whole vector space is re-initialized,
  // which is fine since initialization order for non-essential
  // spaces is irrelevant.
  // Also set _essentialComponent to true for all fully essential spaces.
  std::vector<feSpace *> allEssentialSpaces = _essentialSpaces;
  for(auto *space : _spaces) {
    if(space->getNumComponents() > 1) {
      bool toAdd = false;
      for(int i = 0; i < space->getNumComponents(); ++i) {
        if(space->isEssentialComponent(i)) {
          toAdd = true;
        }
      }
      if(toAdd) allEssentialSpaces.push_back(space);
    }
  }
  for(auto *space : _essentialSpaces) {
    space->setEssentialComponent(0, true);
    space->setEssentialComponent(1, true);
    space->setEssentialComponent(2, true);
  }

  for(feSpace *fS : allEssentialSpaces) {
    if(fS->getDOFInitialization() == dofInitialization::PREVIOUS_SOL) {
      continue;
    }

    int nElm = fS->getNumElements();
    std::vector<feInt> adr(fS->getNumFunctions());
    std::vector<double> coor = fS->getLcoor();
    int numVerticesPerElem = mesh->getCncGeoByName(fS->getCncGeoID())->getNumVerticesPerElem();
    std::vector<double> localCoord(3 * numVerticesPerElem);

    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());

    if(fS->getDOFInitialization() == dofInitialization::NODEWISE ||
       fS->getDOFInitialization() == dofInitialization::EXTRAPOLATED_EULER_0D)
    {
      // Node based: the initial condition is imposed at the vertices DOF
      for(int iElm = 0; iElm < nElm; ++iElm) {
        fS->initializeAddressingVector(iElm, adr);
        mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

        for(int j = 0; j < fS->getNumFunctions(); ++j) {
          double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
          geoSpace->interpolateVectorField(localCoord, r, args.pos);

          if(fS->getNumComponents() > 1) {
            // Vector valued finite element space
            fS->evalFun(args, vecVal);
            int indexComponent = j % fS->getNumComponents();
            val = vecVal[indexComponent];
            if(fS->isEssentialComponent(indexComponent)) {
              _sol[adr[j]] = val;
              if(solContainer != nullptr) {
                solContainer->_sol[0][adr[j]] = val;
              };
            }

            ///////////////////////////////////////////////
            // // Solution time derivative
            // const double x = args.pos[0];
            // const double y = args.pos[1];
            // const double t = args.t;
            // const double a = 0.;
            // const double b = 0.;
            // // vecVal[0] = 5.*t*t*t*t * (a * pow(x, 3) + x * y*y);
            // // vecVal[1] = 5.*t*t*t*t * b * x*y;
            // // vecVal[2] = 0.;

            // // vecVal[0] = -exp(-t) * (a * pow(x, 3) + x * y*y);
            // // vecVal[1] = -exp(-t) * b * x*y;
            // // vecVal[2] = 0.;

            // // vecVal[0] = -2./((1.+t)*(1.+t)) * (a * pow(x, 3) + x * y*y);
            // // vecVal[1] = -2./((1.+t)*(1.+t)) * b * x*y;
            // // vecVal[2] = 0.;

            // vecVal[0] = -2./((1.+t)*(1.+t)) * (a * y + x);
            // vecVal[1] = -2./((1.+t)*(1.+t)) * b * x*y;
            // vecVal[2] = 0.;

            // // const double t = args.t;
            // // vecVal[0] = -exp(-t);
            // // vecVal[1] = -exp(-t);
            // // vecVal[2] = 0.;

            // indexComponent = j % fS->getNumComponents();
            // val = vecVal[indexComponent];
            // if(fS->isEssentialComponent(indexComponent)) {
            //   _dsoldt[adr[j]] = val;
            //   if(solContainer != nullptr) {
            //     solContainer->_solDot[0][adr[j]] = val;
            //   };
            // }
            ///////////////////////////////////////////////

          } else {
            // Scalar valued finite element space
            val = fS->evalFun(args);
            _sol[adr[j]] = val;
            if(solContainer != nullptr) {
              solContainer->_sol[0][adr[j]] = val;
            };

            // ///////////////////////////////////////////////
            // // const double x = args.pos[0];
            // // const double t = args.t;
            // // const double udot = -2./((1.+t)*(1.+t)) * x;
            // // const double udot = -2./((1.+t)*(1.+t));
            // const double udot = (double) rand() / (double)RAND_MAX;
            // // const double udot = nan("");
            // _dsoldt[adr[j]] = udot;
            // if(solContainer != nullptr) {
            //   solContainer->_solDot[0][adr[j]] = udot;
            // };
            // ///////////////////////////////////////////////

          }
        }
      }
    } else if(fS->getDOFInitialization() == dofInitialization::LEAST_SQUARES) {
      // Not yet implemented
      feErrorMsg(FE_STATUS_ERROR, "Cannot initialize essential DOF in the LEAST_SQUARES sense.");
      exit(-1);
    }
  }
}

void feSolution::copySpace(feMesh *mesh, feSpace *s1, feSpace *s2)
{
  int nElm = mesh->getNumElements(s1->getCncGeoID());
  std::vector<feInt> adr1(s1->getNumFunctions());
  std::vector<feInt> adr2(s2->getNumFunctions());
  std::vector<double> coor = s1->getLcoor();
  int numVerticesPerElem = mesh->getCncGeoByName(s1->getCncGeoID())->getNumVerticesPerElem();
  std::vector<double> localCoord(3 * numVerticesPerElem);
  std::vector<double> x(3);

  for(int iElm = 0; iElm < nElm; ++iElm) {
    s1->initializeAddressingVector(iElm, adr1);
    s2->initializeAddressingVector(iElm, adr2);
    for(int j = 0; j < s1->getNumFunctions(); ++j) {
      _sol[adr2[j]] = _sol[adr1[j]];
    }
  }
}

// Adds coeff * the value of sourceSpace into targetSpace.
// Interpolation is required if both spaces do not have the same number of DOFs
// Both spaces should have the same number of components.
feStatus feSolution::addConstantTimesSpace(feMesh *mesh, 
  const double coeff, feSpace *sourceSpace, feSpace *targetSpace)
{
  int nElm = mesh->getNumElements(sourceSpace->getCncGeoID());
  int nCompS = sourceSpace->getNumComponents();
  int nCompT = targetSpace->getNumComponents();

  if(nCompT != nCompS) {
    return feErrorMsg(FE_STATUS_ERROR, "Cannot add space because target space does not have the same number of components!");
  }

  int ndofS = sourceSpace->getNumFunctions();
  int ndofT = targetSpace->getNumFunctions();
  bool interpolate = (ndofS != ndofT);
  std::vector<feInt> adrS(ndofS);
  std::vector<feInt> adrT(ndofT);

  // Very ugly: keep a set of already incremented DOFs...
  std::set<feInt> incrementedDOFs;

  if(interpolate) {
    
    const std::vector<double> &Lcoor = targetSpace->getLcoor();
    double r[3];
    std::vector<double> localSol(ndofS);

    // Add coeff * interpolated field 
    for(int iElm = 0; iElm < nElm; ++iElm) {
      sourceSpace->initializeAddressingVector(iElm, adrS);
      targetSpace->initializeAddressingVector(iElm, adrT);
      this->getSolAtDOF(adrS, localSol);

      // Interpolate each component
      for(int j = 0; j < nCompS; ++j) {
        for(int k = 0; k < ndofT / nCompT; ++k) {
          if(incrementedDOFs.find(adrT[nCompT * k + j]) == incrementedDOFs.end()) {
            r[0] = Lcoor[3 * nCompT * k + 3 * j + 0];
            r[1] = Lcoor[3 * nCompT * k + 3 * j + 1];
            r[2] = Lcoor[3 * nCompT * k + 3 * j + 2];
            double res = sourceSpace->interpolateVectorFieldComponent(localSol, j, r);
            _sol[adrT[nCompT * k + j]] += coeff * res;
            incrementedDOFs.insert(adrT[nCompT * k + j]);
          }
        }
      }
    }

  } else {
    // Add coeff * field at each DOF
    for(int iElm = 0; iElm < nElm; ++iElm) {
      sourceSpace->initializeAddressingVector(iElm, adrS);
      targetSpace->initializeAddressingVector(iElm, adrT);
      for(int j = 0; j < ndofS; ++j) {
        if(incrementedDOFs.find(adrT[j]) == incrementedDOFs.end()) {
          _sol[adrT[j]] += coeff * _sol[adrS[j]];
          incrementedDOFs.insert(adrT[j]);
        }
      }
    }
  }

  return FE_STATUS_OK;
}

// Adds coeff * the squared norm of a space sourceSpace into targetSpace.
// The target space must have as many DOFs as the number of DOFs of a 
// scalar component of the source space.
feStatus feSolution::addSquaredNormOfVectorSpace(feMesh *mesh, 
  const double coeff, feSpace *sourceSpace, feSpace *targetSpace)
{
  int nElm = mesh->getNumElements(sourceSpace->getCncGeoID());
  int nCompS = sourceSpace->getNumComponents();
  int nCompT = targetSpace->getNumComponents();

  if(nCompT != 1) {
    return feErrorMsg(FE_STATUS_ERROR,
      "Source space has %d components and target space has %d components. "
      "Cannot copy norm into target space because target space is not a scalar FE space!\n",
      nCompS, nCompT);
  }

  // Very ugly: keep a set of already incremented DOFs...
  std::set<feInt> incrementedDOFs;

  int ndofS = sourceSpace->getNumFunctions();
  int ndofT = targetSpace->getNumFunctions();
  std::vector<feInt> adrS(ndofS);
  std::vector<feInt> adrT(ndofT);

  for(int iElm = 0; iElm < nElm; ++iElm) {
    sourceSpace->initializeAddressingVector(iElm, adrS);
    targetSpace->initializeAddressingVector(iElm, adrT);
    for(int j = 0; j < ndofT; ++j) {
      double normSquared = 0.;
      for(int k = 0; k < nCompS; ++k) {
        normSquared += _sol[adrS[nCompS * j + k]] * _sol[adrS[nCompS * j + k]];
      }
      if(incrementedDOFs.find(adrT[j]) == incrementedDOFs.end()) {
        _sol[adrT[j]] += coeff * normSquared;
        incrementedDOFs.insert(adrT[j]);
      }
    }
  }

  return FE_STATUS_OK;
}



void feSolution::setSolFromContainer(const feSolutionContainer *solContainer, const int iSol)
{
  const std::vector<double> &solFromContainer = solContainer->getSolution(iSol);
  if(_sol.size() != solFromContainer.size()) {
    _sol.resize(solFromContainer.size());
    _dsoldt.resize(solFromContainer.size());
  }
  size_t nDOFs = solFromContainer.size();
  for(size_t i = 0; i < nDOFs; ++i) {
    _sol[i] = solFromContainer[i];
  }
}

void feSolution::setSolDotToZero()
{
  size_t nDOFs = _dsoldt.size();
  for(size_t i = 0; i < nDOFs; ++i) _dsoldt[i] = 0.0;
}

// Only prints the current solution, not the full container
// The export function should be for the container
// to export the full solution history
void feSolution::printSol(std::string file)
{
  FILE *f = nullptr;
  if(file != "") {
    f = fopen(file.c_str(), "w");
    fprintf(f, "%12.16f\n", _tn); // Print the current time
    fprintf(f, "%ld\n", _sol.size()); // Print the number of DOFs
  }
  for(auto const &val : _sol) {
    if(file != "") {
      fprintf(f, "%12.16f\n", val);
    } else {
      printf("%12.16f\n", val);
    }
  }
  if(file != "") {
    // Print soldot
    for(auto const &val : _dsoldt) {
      fprintf(f, "%12.16f\n", val);
    }
    fclose(f);
  }
}
