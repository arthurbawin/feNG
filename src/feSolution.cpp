#include "feSolution.h"

#include <iostream>
#include <fstream>

feSolution::feSolution(feMesh *mesh, const std::vector<feSpace *> &space,
                       const std::vector<feSpace *> &essBC, feMetaNumber *metaNumber)
  : _dim(mesh->getDim()), _space(space), _essBC(essBC), _c0(0.), _tn(0.)
{
  _sol.resize(metaNumber->getNbDOFs());
  _dsoldt.resize(metaNumber->getNbDOFs());
  for(int i = 0; i < metaNumber->getNbDOFs(); ++i) {
    _sol[i] = 0.0;
    _dsoldt[i] = 0.0;
  }
}

/* Constructs an feSolution from a file created by feSolution::printSol. */
feSolution::feSolution(std::string solutionFile)
  : _dim(0), _space(std::vector<feSpace *>()), _essBC(std::vector<feSpace *>()), _c0(0.)
{
  printf("Info feSolution::feSolution : Reading solution file : %s\n", solutionFile.c_str());
  std::filebuf fb;
  if(fb.open(solutionFile, std::ios::in)) {
    std::istream input(&fb);
    std::string buffer;
    double solutionTime;
    int solutionNDOFs;
    input >> solutionTime;
    getline(input, buffer);
    input >> solutionNDOFs;

    _tn = solutionTime;

    _sol.resize(solutionNDOFs);
    _dsoldt.resize(solutionNDOFs);

    double val;
    for(int i = 0; i < solutionNDOFs; ++i) {
      getline(input, buffer);
      input >> val;
      _sol[i] = val;
    }
    for(int i = 0; i < solutionNDOFs; ++i) {
      getline(input, buffer);
      input >> val;
      _dsoldt[i] = val;
    }

    fb.close();
  } else {
    printf("In feSolution::feSolution : Error - Solution could not be opened.\n");
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
  for(feSpace *const &fS : _space) {
    // std::cout<<"Space "<<fS->getFieldID()<<" - "<<fS->getCncGeoID()<<std::endl;
    int nElm = mesh->getNbElm(fS->getCncGeoID()); // On pourrait donner un _nElm au feSpace
    std::vector<double> coor = fS->getLcoor();
    // std::cout<<"coor"<<std::endl;
    // for(auto const &val : coor)
    // 	std::cout<<val<<std::endl;
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    // std::cout<<geoSpace->getNbFunctions()<<"?"<<std::endl;
    for(int iElm = 0; iElm < nElm; ++iElm) {
      // Call initializeAddressingVector() on each element of each feSpace
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);

      std::vector<double> localCoord = mesh->getCoord(fS->getCncGeoID(), iElm);
      // std::cout<<"size"<<localCoord.size()<<std::endl;

      for(int j = 0; j < fS->getNbFunctions(); ++j) {
        double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
        std::vector<double> x(3, 0.0);
        // x[0] = geoSpace->interpolateField(localCoord, r);
        geoSpace->interpolateVectorField(localCoord, r, x);

        // If the function given to the feSpace is null, the solution is not changed,
        // i.e. we continue with the value from a previous computation.
        if(fS->isFctDefined()) {
          double val = fS->evalFun(_tn, x);
          _sol[fS->getAddressingVectorAt(j)] = val;
        } else {
          // printf("Sol un changed in field %s - %s\n", fS->getFieldID().c_str(),
          // fS->getCncGeoID().c_str());
        }
      }
    }
  }
}

void feSolution::initializeEssentialBC(feMesh *mesh, feMetaNumber *metaNumber,
                                       feSolutionContainer *solContainer)
{
  for(feSpace *const &fS : _essBC) {
    int nElm = mesh->getNbElm(fS->getCncGeoID());
    // std::cout<<"Initializing essential BC on "<<nElm<<" elements on cnc
    // "<<fS->getCncGeoID()<<std::endl;
    std::vector<double> coor = fS->getLcoor();
    feSpace *geoSpace = mesh->getGeometricSpace(fS->getCncGeoID());
    for(int iElm = 0; iElm < nElm; ++iElm) {
      fS->initializeAddressingVector(metaNumber->getNumbering(fS->getFieldID()), iElm);
      // std::cout<<"On elem "<<iElm<<std::endl;
      for(int kk = 0; kk < fS->getNbFunctions(); ++kk) {
        // std::cout<<fS->getAddressingVectorAt(kk)<<std::endl;
      }
      std::vector<double> localCoord = mesh->getCoord(fS->getCncGeoID(), iElm);
      for(int j = 0; j < fS->getNbFunctions(); ++j) {
        double r[3] = {coor[3 * j], coor[3 * j + 1], coor[3 * j + 2]};
        std::vector<double> x(3, 0.0);
        geoSpace->interpolateVectorField(localCoord, r, x);
        if(fS->isFctDefined()) {
          _sol[fS->getAddressingVectorAt(j)] = fS->evalFun(_tn, x);
          // printf("Initializing _sol[%d] = %f en (%+-4.4e, %+-4.4e)\n",
          // fS->getAddressingVectorAt(j), _sol[fS->getAddressingVectorAt(j)], x[0], x[1]);
        } else {
          printf("BC Sol un changed in field %s - %s\n", fS->getFieldID().c_str(),
                 fS->getCncGeoID().c_str());
        }
        // printf("_sol[%d] = %f\n", fS->getAddressingVectorAt(j),
        // _sol[fS->getAddressingVectorAt(j)]);
        if(solContainer != nullptr) {
          if(fS->isFctDefined()) {
            solContainer->_sol[0][fS->getAddressingVectorAt(j)] = fS->evalFun(_tn, x);
          } else {
            // printf("BC Sol un changed in field %s - %s\n", fS->getFieldID().c_str(),
            // fS->getCncGeoID().c_str());
          }
        };
      }
    }
  }
}

void feSolution::setSolFromContainer(feSolutionContainer *solContainer, int iSol)
{
  std::vector<double> solFromContainer = solContainer->getSolution(iSol);
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
    fprintf(f, "%d\n", _sol.size()); // Print the number of DOFs
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
