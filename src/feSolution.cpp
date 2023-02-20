#include "feSolution.h"

#include <iostream>
#include <fstream>

feSolution::feSolution(int numDOF, const std::vector<feSpace *> &space,
                       const std::vector<feSpace *> &essentialSpaces)
  : _spaces(space), _essentialSpaces(essentialSpaces), _nDOF(numDOF), _c0(0.), _tn(0.)
{
  _sol.resize(numDOF);
  _dsoldt.resize(numDOF);
  for(int i = 0; i < numDOF; ++i) {
    _sol[i] = 0.0;
    _dsoldt[i] = 0.0;
  }
}

/* Constructs an feSolution from a file created by feSolution::printSol. */
feSolution::feSolution(std::string solutionFile)
  : _spaces(std::vector<feSpace *>()), _essentialSpaces(std::vector<feSpace *>()), _c0(0.)
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

void feSolution::initializeTemporalSolution(double t0, double t1, int nTimeSteps)
{
  _t0 = t0;
  _t1 = t1;
  _nTimeSteps = nTimeSteps;
  _dt = (t1 - t0) / (double)nTimeSteps;
  _tn = t0;
}

void feSolution::initializeUnknowns(feMesh *mesh)
{
  std::vector<double> x(3);
  std::vector<double> vecVal(3);
  double val;

  for(feSpace *fS : _spaces) {
    if(fS->getDOFInitialization() == dofInitialization::PREVIOUS_SOL) {
      continue;
    }

    int nElm = fS->getNumElements();
    std::vector<feInt> adr(fS->getNumFunctions());
    std::vector<double> coor = fS->getLcoor();
    int numVerticesPerElem = mesh->getCncGeoByName(fS->getCncGeoID())->getNumVerticesPerElem();
    std::vector<double> localCoord(3 * numVerticesPerElem);

    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());

    if(fS->getDOFInitialization() == dofInitialization::NODEWISE) {
      // Node based: the initial condition is imposed at the vertices DOF
      for(int iElm = 0; iElm < nElm; ++iElm) {
        fS->initializeAddressingVector(iElm, adr);
        mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

        for(int j = 0; j < fS->getNumFunctions(); ++j) {
          double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
          geoSpace->interpolateVectorField(localCoord, r, x);

          if(fS->getNumComponents() > 1) {
            // Vector valued finite element space
            fS->evalFun(_tn, x, vecVal);
            int indexComponent = j % fS->getNumComponents();
            val = vecVal[indexComponent];
          } else {
            // Scalar valued finite element space
            val = fS->evalFun(_tn, x);
          }

          _sol[adr[j]] = val;
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
            geoSpace->interpolateVectorField(localCoord, r, x);
            uQuad = fS->evalFun(_tn, x);
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
  }
}

void feSolution::initializeEssentialBC(feMesh *mesh, feSolutionContainer *solContainer)
{
  std::vector<double> x(3);
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
       fS->getDOFInitialization() == dofInitialization::EXTRAPOLATED_EULER_0D) {
      // Node based: the initial condition is imposed at the vertices DOF
      for(int iElm = 0; iElm < nElm; ++iElm) {
        fS->initializeAddressingVector(iElm, adr);
        mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

        for(int j = 0; j < fS->getNumFunctions(); ++j) {
          double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
          geoSpace->interpolateVectorField(localCoord, r, x);

          if(fS->getNumComponents() > 1) {
            // Vector valued finite element space
            fS->evalFun(_tn, x, vecVal);
            int indexComponent = j % fS->getNumComponents();
            val = vecVal[indexComponent];
            if(fS->isEssentialComponent(indexComponent)) {
              _sol[adr[j]] = val;
              if(solContainer != nullptr) {
                solContainer->_sol[0][adr[j]] = val;
              };
            }
          } else {
            // Scalar valued finite element space
            val = fS->evalFun(_tn, x);
            _sol[adr[j]] = val;
            if(solContainer != nullptr) {
              solContainer->_sol[0][adr[j]] = val;
            };
          }
        }
      }
    } else if(fS->getDOFInitialization() == dofInitialization::LEAST_SQUARES) {
      // Not yet implemented
      feErrorMsg(FE_STATUS_ERROR, "Cannot initialize essential DOF in the LEAST_SQUARES sense.");
      exit(-1);
    }
  }

  // if(fS->getDOFInitialization() == dofInitialization::EXTRAPOLATED_EULER_0D) {
  //   // Extrapolate solution from inside the domain to the boundary
  //   // Only in 1D for Euler equations for now

  //   for(int iElm = 0; iElm < nElm; ++iElm) {
  //     fS->initializeAddressingVector(iElm, adr);

  //     for(int j = 0; j < fS->getNumFunctions(); ++j) {
  //       double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
  //       geoSpace->interpolateVectorField(localCoord, r, x);

  //       // Extrapolate solution

  //       _sol[adr[j]] = val;
  //     }

  //   }
  // }
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

void feSolution::setSolFromContainer(feSolutionContainer *solContainer, int iSol)
{
  std::vector<double> &solFromContainer = solContainer->getSolution(iSol);
  if(_sol.size() != solFromContainer.size()) {
    _sol.resize(solFromContainer.size());
    _dsoldt.resize(solFromContainer.size());
  }
  size_t nDOFs = solFromContainer.size();
  for(size_t i = 0; i < nDOFs; ++i) _sol[i] = solFromContainer[i];
}

void feSolution::setSolDotToZero()
{
  size_t nDOFs = _dsoldt.size();
  for(size_t i = 0; i < nDOFs; ++i) _dsoldt[i] = 0.0;
}

void feSolution::printSol(std::string file)
{
  FILE *f;
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
