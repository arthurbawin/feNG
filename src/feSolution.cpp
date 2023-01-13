#include "feSolution.h"

#include <iostream>
#include <fstream>

feSolution::feSolution(int numDOF, const std::vector<feSpace *> &space,
                       const std::vector<feSpace *> &essentialSpaces)
  : _spaces(space)
  , _essentialSpaces(essentialSpaces)
  , _nDOF(numDOF)
  , _c0(0.)
  , _tn(0.)
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
  : _spaces(std::vector<feSpace *>()), 
  _essentialSpaces(std::vector<feSpace *>()),
  _c0(0.)
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

void feSolution::initializeTemporalSolution(double t0, double t1, int nTimeSteps)
{
  _t0 = t0;
  _t1 = t1;
  _nTimeSteps = nTimeSteps;
  _dt = (t1 - t0) / (double)nTimeSteps;
  _tn = t0;
}

void feSolution::initializeUnknowns(feMesh *mesh, feMetaNumber *metaNumber)
{
  for(feSpace *const &fS : _spaces) {

    int nElm = mesh->getNbElm(fS->getCncGeoID()); // On pourrait donner un _nElm au feSpace
    std::vector<feInt> adr(fS->getNbFunctions());
    std::vector<double> coor = fS->getLcoor();
    int nbNodePerElem = mesh->getCncGeoByName(fS->getCncGeoID())->getNbNodePerElem();
    std::vector<double> localCoord(3*nbNodePerElem);
    std::vector<double> x(3);
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    
    if(true){
      // Node based: the initial condition is imposed at the vertices DOF
      for(int iElm = 0; iElm < nElm; ++iElm) {
        fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm, adr);
        mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);
        
        for(int j = 0; j < fS->getNbFunctions(); ++j) {
          double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
          // x[0] = geoSpace->interpolateField(localCoord, r);
          geoSpace->interpolateVectorField(localCoord, r, x);

          // If the function given to the feSpace is null, the solution is not changed,
          // i.e. we continue with the value from a previous computation.
          if(fS->isFctDefined()) {
            double val = fS->evalFun(_tn, x);
            _sol[adr[j]] = val;
          }
        }
      }
    } else{
      // Cell based 1D Legendre polynomials: the initial condition is satisfied on each element in the least-square sense:
      //
      //  <phi_i, phi_j> U_j^e = <phi_i, funSol> ==> A_ij^e U_j^e = B_j^e
      //
      // !!! We assume that the jacobian cancels out from both side, which is true only for P1 line geometry !!!

      // The analytic diagonal mass matrix A_ii = 2/(2*i+1)
      std::vector<double> A(fS->getNbFunctions(), 0.);
      std::vector<double> B(fS->getNbFunctions(), 0.);
      for(int i = 0; i < fS->getNbFunctions(); ++i){
        A[i] = 2./(2. * i + 1.);
      }

      int nQuad = geoSpace->getNbQuadPoints();
      std::vector<double> wQuad = geoSpace->getQuadratureWeights();
      std::vector<double> xQuad = geoSpace->getRQuadraturePoints();

      // std::vector<double> &J = mesh->getCncGeoByName(fS->getCncGeoID())->getJacobians();

      double uQuad = 0.;
      std::vector<double> phi(fS->getNbFunctions(), 0.);
      for(int iElm = 0; iElm < nElm; ++iElm){

        fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm, adr);
        mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);

        // Compute the RHS
        for(int i = 0; i < fS->getNbFunctions(); ++i){

          B[i] = 0.;

          for(int k = 0; k < nQuad; ++k){

            // Evaluate analytic solution at quad points
            double r[3] = {xQuad[k], 0., 0.};
            geoSpace->interpolateVectorField(localCoord, r, x);
            uQuad = fS->evalFun(_tn, x);
            // Evaluate Legendre polynomials at quad points
            fS->L(r, phi.data());

            // feInfo("u = %f - phi(%f) = (%f - %f - %f vs %f)", uQuad, r[0], phi[0], phi[1], phi[2], 3./2.*r[0]*r[0] - 0.5);

            B[i] += wQuad[k] * phi[i] * uQuad;
          }
        }

        // Solve Aij Uj = Bj on the element
        std::vector<double> elmSolution(adr.size());
        for(int i = 0; i < fS->getNbFunctions(); ++i){
          _sol[adr[i]] = B[i] / A[i];
          elmSolution[i] = _sol[adr[i]];
          // feInfo("elm %d - %d - dof %d: Bi = %f Ai = %f - Computed initial condition %f", iElm, i, adr[i], B[i], A[i], B[i] / A[i];
        }

        double rr[3] = {-1., 0., 0.};
        geoSpace->interpolateVectorField(localCoord, rr, x);
        double myVal = fS->interpolateField(elmSolution, rr);
        feInfo("Interpolation en x = %f : uh = %f - u = %f - adr = %d %d %d", x[0], myVal, fS->evalFun(_tn, x), adr[0], adr[1], adr[2]);

      }
    }
  }
}

void feSolution::initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber,
                                       feSolutionContainer *solContainer)
{
  for(feSpace *const &fS : _essentialSpaces) {
    std::vector<feInt> adrS(fS->getNbFunctions());
    int nElm = mesh->getNbElm(fS->getCncGeoID());
    std::vector<double> coor = fS->getLcoor();
    int nbNodePerElem = mesh->getCncGeoByName(fS->getCncGeoID())->getNbNodePerElem();
    std::vector<double> localCoord(3*nbNodePerElem); 
    std::vector<double> x(3);
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    
    for(int iElm = 0; iElm < nElm; ++iElm) {
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm, adrS);
      mesh->getCoord(fS->getCncGeoID(), iElm, localCoord);
      for(int j = 0; j < fS->getNbFunctions(); ++j) {
        double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
        geoSpace->interpolateVectorField(localCoord, r, x);
        
        if(fS->isFctDefined()) {
          _sol[adrS[j]] = fS->evalFun(_tn, x);
          // fS->getAddressingVectorAt(j), _sol[fS->getAddressingVectorAt(j)], x[0], x[1]);
        } else {
          printf("BC Sol un changed in field %s - %s\n", fS->getFieldID().c_str(),
                 fS->getCncGeoID().c_str());
        }
        // printf("_sol[%d] = %f\n", fS->getAddressingVectorAt(j),
        // _sol[fS->getAddressingVectorAt(j)]);
        if(solContainer != nullptr) {
          if(fS->isFctDefined()) {
            solContainer->_sol[0][adrS[j]] = fS->evalFun(_tn, x);
          } else {
            // printf("BC Sol un changed in field %s - %s\n", fS->getFieldID().c_str(),
            // fS->getCncGeoID().c_str());
          }
        };
      }
    }
  }
}

void feSolution::copySpace(feMesh *mesh, feMetaNumber *metaNumber, feSpace *s1, feSpace *s2)
{
  int nElm = mesh->getNbElm(s1->getCncGeoID());
  std::vector<feInt> adr1(s1->getNbFunctions());
  std::vector<feInt> adr2(s2->getNbFunctions());
  std::vector<double> coor = s1->getLcoor();
  int nbNodePerElem = mesh->getCncGeoByName(s1->getCncGeoID())->getNbNodePerElem();
  std::vector<double> localCoord(3*nbNodePerElem);
  std::vector<double> x(3);
  
  for(int iElm = 0; iElm < nElm; ++iElm) {
    s1->initializeAddressingVector(metaNumber->getNumbering(s1->getFieldID()), iElm, adr1);
    s2->initializeAddressingVector(metaNumber->getNumbering(s2->getFieldID()), iElm, adr2);    
    for(int j = 0; j < s1->getNbFunctions(); ++j) {
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
