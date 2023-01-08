#include "feRecovery.h"
#include "feNG.h"

#include <algorithm>
#include "../contrib/Eigen/Dense"

static bool isBoundary(feMesh *mesh, int vertex)
{
  // return false;
  return fabs(mesh->getVertex(vertex)->x()) < 1e-10 ||
         fabs(mesh->getVertex(vertex)->x() - 5.0) < 1e-10;
}

static bool isBoundary(feMesh *mesh, Edge edge) { return false; }

static double boundaryCondition(feMesh *mesh, int vertex, feFunction *solRef)
{
  return solRef->eval(
    0, {mesh->getVertex(vertex)->x(), mesh->getVertex(vertex)->y(), mesh->getVertex(vertex)->z()});
}

static void boundaryConditionVec(feMesh *mesh, int vertex, feVectorFunction *solRefGrad,
                                 std::vector<double> &res)
{
  solRefGrad->eval(
    0, {mesh->getVertex(vertex)->x(), mesh->getVertex(vertex)->y(), mesh->getVertex(vertex)->z()},
    res);
}

static double boundaryCondition(feMesh *mesh, Edge edge, feFunction *solRef) { return 0.0; }

static void boundaryConditionVec(feMesh *mesh, Edge edge, feVectorFunction *solRefGrad,
                                 std::vector<double> &res)
{
}

static inline double matNorm2(const std::vector<double> &v1, const std::vector<double> &v2, int n)
{
  double sqr = 0;
  for(int i = 0; i < n; i++) {
    sqr += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  }
  return sqrt(sqr);
}

static inline double myPow(double base, int exp){
  // if (exp == 0)
  //     return 1.;
  // else if (exp % 2)
  //     return base * myPow(base, exp - 1);
  // else {
  //     int temp = myPow(base, exp / 2);
  //     return temp * temp;
  // }

  // return pow(base, exp);
  double res = 1.;
  for(int i = 0; i < exp; ++i)
    res *= base;
  return res;
}

// The names of the reconstructed fields
static std::map<std::pair<int, int>, std::string> suffix = {
  {{0, 0}, ""},     {{1, 0}, "dx"},   {{1, 1}, "dy"},   {{2, 0}, "dxx"},  {{2, 1}, "dxy"},
  {{2, 2}, "dyx"},  {{2, 3}, "dyy"},  {{3, 0}, "dxxx"}, {{3, 1}, "dxxy"}, {{3, 2}, "dxyx"},
  {{3, 3}, "dxyy"}, {{3, 4}, "dyxx"}, {{3, 5}, "dyxy"}, {{3, 6}, "dyyx"}, {{3, 7}, "dyyy"},
};

fePatch::fePatch(feCncGeo *cnc, feMesh *mesh)
{
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getNodeConnectivityCopy();
  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  _nVertices = _vertices.size();
  _nNodePerElm = cnc->getNbNodePerElem();
  _nEdgePerElm = cnc->getNbEdgePerElem();
  std::vector<int> &connecNodes = cnc->getNodeConnectivityRef();

  int nElm = cnc->getNbElm();
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < _nNodePerElm; ++j) {
      vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i); // TODO : A verifier
    }
  }

  switch(mesh->getDim()) {
    case 1:
      // The 1D patch associated to vertices is the 2 surrounding elements.
      // The patch associated to edges is the union of the extremities' patches.
      for(auto e : mesh->_edges) {
        edgeToElems[e.getTag()].insert(vertToElems[e.getTag(0)].begin(),
                                       vertToElems[e.getTag(0)].end());
        edgeToElems[e.getTag()].insert(vertToElems[e.getTag(1)].begin(),
                                       vertToElems[e.getTag(1)].end());
      }
      // Increase the size of the patch if there are too few elements
      for(auto &v : vertToElems) {
        if(v.second.size() < 2) {
          std::cout << "v " << v.first << std::endl;
          std::cout << std::endl;
          for(auto &val : v.second) std::cout << "elem " << val << std::endl;
          std::cout << std::endl;
          std::vector<int> toAdd;
          for(auto val : toAdd) std::cout << "toAdd " << val << std::endl;
          std::cout << std::endl;
          for(auto e : v.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto val : toAdd) std::cout << "toAdd " << val << std::endl;
          std::cout << std::endl;
          for(auto e : toAdd) v.second.insert(e);
        }
      }
      for(auto &v : vertToElems) {
        std::cout << "v " << v.first << std::endl;
        for(auto &val : v.second) std::cout << "elem " << val << std::endl;
        std::cout << std::endl;
      }

      break;
    case 2:

      for(auto e : mesh->_edges) {
        // Insert the patches of both edge vertices to the edge's patch
        int v0 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(0));
        int v1 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(1));
        edgeToElems[e.getTag()].insert(vertToElems[v0].begin(), vertToElems[v0].end());
        edgeToElems[e.getTag()].insert(vertToElems[v1].begin(), vertToElems[v1].end());
      }

      // Increase the size of the patch if there are too few elements
      for(auto &v : vertToElems) {
        if(v.second.size() <= 2) {
          std::vector<int> toAdd;
          for(auto e : v.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto e : toAdd) v.second.insert(e);
        }
      }

      for(auto &p : edgeToElems) {
        if(p.second.size() <= 2) {
          std::vector<int> toAdd;
          for(auto e : p.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto e : toAdd) p.second.insert(e);
        }
      }
      break;
    default:
      printf("In fePatch : Error - Element patches are only defined for 1D and 2D meshes.\n");
      return;
  }
}

feRecovery::feRecovery(feMetaNumber *metaNumber, feSpace *space, feMesh *mesh, feSolution *sol,
                       std::vector<double> &norm, feFunction *solRef, std::string meshName,
                       std::string metricMeshName, feVectorFunction *solRefGrad,
                       feVectorFunction *solRefHess, feFunction *fund3udx, bool append)
  : _metaNumber(metaNumber), _mesh(mesh), _sol(sol), _intSpace(space), _solRef(solRef),
    _solRefGrad(solRefGrad), _solRefHess(solRefHess)
{
  _cnc = space->getCncGeo();
  _nElm = _cnc->getNbElm();
  _nNodePerElm = _cnc->getNbNodePerElem();
  _adr.resize(_intSpace->getNbFunctions());
  _solution.resize(_adr.size());
  _geoSpace = _cnc->getFeSpace();
  _degSol = space->getPolynomialDegree();
  _patch = new fePatch(_cnc, _mesh);

  _dim = mesh->getDim();

  // Set total number of recoveries/derivations
  _nTotalRecoveries = pow(_dim, _degSol); // TODO : A verifier
  _nTotalDerivations = pow(_dim, _degSol + 1);

  // The dimension of the polynomial bases for recoveries is the one of degree k+1 :
  if(_dim == 1) {
    _dimRecovery = _degSol + 2;
    _dimDerivation = _degSol + 1;
    _dim2Derivation = _degSol;
  } else if(_dim == 2) {
    if(_cnc->getForme() == "TriP1" || _cnc->getForme() == "TriP2") {
      _dimRecovery = (_degSol + 2) * (_degSol + 3) / 2;
      _dimDerivation = (_degSol + 1) * (_degSol + 2) / 2;
    } else if(_cnc->getForme() == "QuadP1" || _cnc->getForme() == "QuadP2") {
      _dimRecovery = (_degSol + 2) * (_degSol + 2);
      _dimDerivation = (_degSol + 1) * (_degSol + 1);
    } else {
      printf("Error : Mesh connectivity \"%s\" not available for polynomial recovery.\n",
             _cnc->getForme().c_str());
    }
  } else if(_dim == 3) {
    if(_cnc->getForme() == "TetP1" || _cnc->getForme() == "TetP2") {
      _dimRecovery = (_degSol + 2) * (_degSol + 3) * (_degSol + 4) / 6;
      _dimDerivation = (_degSol + 1) * (_degSol + 2) * (_degSol + 3) / 6;
    } else {
      printf("Error : Mesh connectivity \"%s\" not available for polynomial recovery.\n",
             _cnc->getForme().c_str());
    }
  }

  // // Resize coefficients vectors
  // std::vector<int> &vertices = _patch->getVertices();

  // int nVertices = vertices.size();
  // int nEdges = _mesh->_edges.size();

  // // int localEdge = fabs(_cnc->getEdgeConnectivity(elem, iEdge));

  // feInfo("size1 = %d", nVertices * _nTotalRecoveries * _dimRecovery);
  // feInfo("size1 = %d", nVertices * _nTotalDerivations * _dimDerivation);
  // feInfo("size2 = %d", nEdges * _nTotalRecoveries * _dimRecovery);
  // feInfo("size2 = %d", nEdges * _nTotalDerivations * _dimDerivation);

  // recoveryCoeffAtVertices.resize(nVertices * _nTotalRecoveries * _dimRecovery, 0.);
  // derivativeCoeffAtVertices.resize(nVertices * _nTotalDerivations * _dimDerivation, 0.);

  // recoveryCoeffOnEdges2.resize(nEdges * _nTotalRecoveries * _dimRecovery, 0.);
  // derivativeCoeffOnEdges2.resize(nEdges * _nTotalDerivations * _dimDerivation, 0.);

  // The exponents of the monomials :
  _expXRecovery.resize(_dimRecovery, 0);
  _expYRecovery.resize(_dimRecovery, 0);
  _expZRecovery.resize(_dimRecovery, 0);

  _expX.resize(_dimDerivation, 0);
  _expY.resize(_dimDerivation, 0);
  _expZ.resize(_dimDerivation, 0);

  int ind = 0, n = _degSol + 1;

  if(_dim == 1) {
    _expX2Derivation.resize(_dim2Derivation, 0);
    // The 1D basis is 1 x x^2 x^3 ...
    for(int i = 0; i < _dimRecovery; ++i) {
      _expXRecovery[i] = i;
      if(i < _dimDerivation) _expX[i] = i;
      if(i < _dim2Derivation) _expX2Derivation[i] = i;
    }

    // for(auto val : _expXRecovery)
    //   std::cout<<val<<std::endl;

    // std::cout<<std::endl;

    // for(auto val : _expX)
    //   std::cout<<val<<std::endl;
    // std::cout<<std::endl;

  } else if(_dim == 2 && (_cnc->getForme() == "TriP1" || _cnc->getForme() == "TriP2")) {
    for(int j = 0; j <= n; ++j) {
      for(int i = 0; i <= n - j; ++i) {
        _expXRecovery[ind] = i;
        _expYRecovery[ind] = j;
        ++ind;
      }
    }
    ind = 0;
    for(int j = 0; j < n; ++j) {
      for(int i = 0; i < n - j; ++i) {
        _expX[ind] = i;
        _expY[ind] = j;
        ++ind;
      }
    }

    for(auto val : _expXRecovery) std::cout << val << std::endl;

    std::cout << std::endl;

    for(auto val : _expX) std::cout << val << std::endl;
    std::cout << std::endl;

  } else {
    printf(
      "TODO : Error : Recovery coefficients not implemented for chosen dimension and/or element\n");
  }

  // allocateStructures();
  printf("Computing inverses with Eigen...\n");
  switch(_dim) {
    case 1:
      matrixInverseEigen1D();
      break;
    case 2:
      matrixInverseEigen2D();
      break;
  }
  printf("Done\n");

  // To export all derivatives to a file, so that we don't need to store them all
  // std::string derivativesFileName = "derivatives.msh";
  // FILE *dFile = fopen(derivativesFileName.c_str(), "w");

  std::filebuf fbIn, fbOut;
  fbIn.open(meshName, std::ios::in);
  std::istream input(&fbIn);
  fbOut.open(metricMeshName, std::ios::out);
  std::ostream output(&fbOut);

  if(_dim > 1) {
    printf("Info in feRecovery : Derivatives will be written to file \"%s\"\n",
           metricMeshName.c_str());
    std::string buffer;
    if(append) {
      // Copy .msh file
      while(getline(input, buffer)) {
        output << buffer << std::endl;
      }
      fbIn.close();
    } else {
      // Copy .msh file except for the possible previous NodeData
      while(getline(input, buffer)) {
        if(buffer == "$NodeData") {
          while(buffer != "$EndNodeData") getline(input, buffer);
          getline(input, buffer);
        }
        output << buffer << std::endl;
      }
      fbIn.close();
    }
  }

  for(int iDerivative = 0; iDerivative < _degSol + 1; ++iDerivative) {
    // for(int iDerivative = 0; iDerivative < 1; ++iDerivative) {
    // bool recoverDerivative = (iDerivative > 0);

    // if(iDerivative < 2){
    // Recovery of the solution if iDerivative = 0, of the derivatives if > 0
    switch(_dim) {
      case 1:
        solveLeastSquareEigen1D(pow(_dim, iDerivative), iDerivative);
        break;
      case 2:
        solveLeastSquareEigen2D(pow(_dim, iDerivative), iDerivative);
        break;
    }
    // }

    // if(iDerivative < 2){
    for(int i = 0; i < pow(_dim, iDerivative); ++i) {
      derivative(i, iDerivative, output);
      printf("Computed derivatives of solution of order %d\n", iDerivative);
    }
    // } else{
    //   for(int i = 0; i < pow(_dim, iDerivative); ++i) {
    //     secondDerivative(i, iDerivative, output);
    //   }
    // }

    if(iDerivative == 0 && solRef != nullptr) {
      estimateError(norm, solRef);
    }
    if(iDerivative == 0 && solRefGrad != nullptr) {
      estimateH1Error(norm, solRefGrad);
    }
    if(iDerivative == 1 && solRefGrad != nullptr) {
      estimateDudxError(norm, solRefGrad);
    }
    if(iDerivative == 1 && solRefHess != nullptr) {
      estimateHessError(norm, solRefHess);
    }
    if(iDerivative == 2 && fund3udx != nullptr) {
      estimated3Error(norm, fund3udx);
    }
  }

  if(_dim > 1) fbOut.close();

  getErrorPolynomials();
}

// Create a feRecovery from a file in which there are the derivatives at the DOFs
feRecovery::feRecovery(feSpace *space, feMesh *mesh, std::string recoveryFile)
  : _mesh(mesh), _intSpace(space)
{
  _cnc = space->getCncGeo();
  _nElm = _cnc->getNbElm();
  _nNodePerElm = _cnc->getNbNodePerElem();
  _geoSpace = _cnc->getFeSpace();
  _degSol = space->getPolynomialDegree();
  _dim = mesh->getDim();
  _patch = new fePatch(_cnc, _mesh);

  // Read the recovery file
  std::filebuf fb;
  feInfo("Reading recovery file : %s", recoveryFile.c_str());
  std::ifstream f(recoveryFile.c_str());
  if(!f.good()) {
    feWarning("Recovery file does not exist. Cannot create feRecovery.");
  }

  int nRecoveries, nDOFs;
  if(fb.open(recoveryFile, std::ios::in)) {
    std::istream input(&fb);
    std::string buffer;

    // getline(input, buffer);
    input >> nRecoveries;
    feInfo("Reading %d rec", nRecoveries);
    derivAtVertices.resize(nRecoveries);
    for(int i = 0; i < nRecoveries; ++i){
      getline(input, buffer);
      input >> nDOFs;
      feInfo("Reading %d DOFs", nDOFs);
      derivAtVertices[i].resize(nDOFs);
      for(int j = 0; j < nDOFs; ++j){
        getline(input, buffer);
        input >> derivAtVertices[i][j];
      }
    }

    getline(input, buffer);
    input >> nRecoveries;
    feInfo("Reading %d rec", nRecoveries);
    derivAtEdges.resize(nRecoveries);
    for(int i = 0; i < nRecoveries; ++i){
      getline(input, buffer);
      input >> nDOFs;
      feInfo("Reading %d DOFs", nDOFs);
      derivAtEdges[i].resize(nDOFs);
      for(int j = 0; j < nDOFs; ++j){
        getline(input, buffer);
        input >> derivAtEdges[i][j];
      }
    }

    fb.close();
  } // if fb.open

}

void feRecovery::allocateStructures()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr;
  // Create Petsc matrix and vectors for least square resolution
  // std::vector<PetscInt> nnz(_dimRecovery, _dimRecovery);
  // ierr = MatCreate(PETSC_COMM_WORLD, &A);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, _dimRecovery, _dimRecovery);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetFromOptions(A);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSeqAIJSetPreallocation(A, 0, nnz.data());
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetUp(A);

  ierr = MatCreate(PETSC_COMM_WORLD, &A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetType(A, MATDENSE);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, _dimRecovery, _dimRecovery);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetUp(A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecCreate(PETSC_COMM_WORLD, &b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetSizes(b, PETSC_DECIDE, _dimRecovery);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(b, &c);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(b, &res);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecZeroEntries(b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set solver and preconditioner
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetType(ksp, KSPGMRES);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPGetPC(ksp, &preconditioner);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = PCSetType(preconditioner, PCILU);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  PetscReal rel_tol = 1e-6;
  PetscReal abs_tol = 1e-12;
  PetscReal div_tol = 1e6;
  PetscInt max_iter = 500;

  ierr = KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol, max_iter);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Trivial indices for Petsc block assignment
  indI = (PetscInt *)malloc(sizeof(PetscInt) * _dimRecovery);
  indJ = (PetscInt *)malloc(sizeof(PetscInt) * _dimRecovery);
  valA = (PetscScalar *)malloc(sizeof(PetscScalar) * _dimRecovery * _dimRecovery);
  valb = (PetscScalar *)malloc(sizeof(PetscScalar) * _dimRecovery);
  for(int i = 0; i < _dimRecovery; ++i) {
    indI[i] = indJ[i] = i;
  }
#endif
}

// void feRecovery::freeStructures() {
//   // Free structures
//   free(valb);
//   free(valA);
//   free(indJ);
//   free(indI);
//   PetscErrorCode ierr;
//   ierr = KSPDestroy(&ksp);
//   CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = VecDestroy(&res);
//   CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = VecDestroy(&c);
//   CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = VecDestroy(&b);
//   CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = MatDestroy(&A);
//   CHKERRABORT(PETSC_COMM_WORLD, ierr);
// }

void feRecovery::matrixInverseEigen1D()
{
  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord(9, 0.), x(3, 0.0), monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  std::cout << nQuad << std::endl;

  std::vector<double> &solVec = _sol->getSolutionReference();

  // Matrices defined on the vertices
  for(auto v : _patch->getVertices()) {
    double xv = _mesh->getVertex(v)->x();
    Eigen::Matrix<double, 3, 3> myMat3 = Eigen::Matrix<double, 3, 3>::Zero();
    Eigen::Matrix<double, 4, 4> myMat4 = Eigen::Matrix<double, 4, 4>::Zero();

    Eigen::Matrix<double, 5, 5> myMatContrainte = Eigen::Matrix<double, 5, 5>::Zero();

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }

      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        double xLoc = x[0] - xv;
        // printf("Vertex %d (%2.2f)- Elem %d - Quad node %d (%2.2f) - xLoc = %2.2f\n",
        //   v, xv, elem, k, x[0], xLoc);

        for(int i = 0; i < _dimRecovery; ++i) {
          monomials[i] = pow(xLoc, _expXRecovery[i]);
          // monomials[i] = myPow(xLoc, _expXRecovery[i]);
        }

        // std::cout<<"Monomes : "<<std::endl;
        // for(auto val : monomials)
        //   std::cout<<val<<std::endl;
        // std::cout<<std::endl;

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol) {
              case 1:
                myMat3(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
              case 2:
                myMat4(i, j) += jac * w[k] * monomials[i] * monomials[j];
                myMatContrainte(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
            }
          }
          // myMatContrainte(_dimRecovery,1) += jac * w[k] * monomials[i];
          // myMatContrainte(i,_dimRecovery) += jac * w[k] * monomials[i];
        }
      }
      // if(isBoundary(v)){
      // Not an integral :
      myMatContrainte(_dimRecovery, 0) = 1.0;
      myMatContrainte(0, _dimRecovery) = 1.0;
      // }
    }
    switch(_degSol) {
      case 1:
        lsInvAtVertices3[v] = myMat3.inverse();
        break;
      case 2:
        lsInvAtVertices4[v] = myMat4.inverse();
        lsInvAtVerticesCont[v] = myMatContrainte.inverse();
        std::cout << "Vertex:" << std::endl;
        std::cout << myMat4 << std::endl;
        std::cout << std::endl;
        std::cout << myMatContrainte << std::endl;
        std::cout << std::endl;
        std::cout << lsInvAtVerticesCont[v] << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        break;
    }
  }
  printf("Least square matrices at vertices : done\n");

  // Matrices defined on the edges
  for(auto e : _mesh->_edges) {
    // TODO : boucler sur le nombre de DOFS par edge, ici on suppose juste un P2 avec 1 dof

    Eigen::Matrix<double, 3, 3> myMat3 = Eigen::Matrix<double, 3, 3>::Zero();
    Eigen::Matrix<double, 4, 4> myMat4 = Eigen::Matrix<double, 4, 4>::Zero();

    Eigen::Matrix<double, 5, 5> myMatContrainte = Eigen::Matrix<double, 5, 5>::Zero();

    std::cout << "==========================" << std::endl;

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }
      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        double xLoc =
          x[0] - (_mesh->getVertex(e.getTag(0))->x() + _mesh->getVertex(e.getTag(1))->x()) / 2.0;
        // printf("Edge %d (%2.2f)- Elem %d - Quad node %d (%2.2f) - xLoc = %2.2f\n",
        //   e.getTag(), (_mesh->getVertex(e.getTag(0))->x() +
        //   _mesh->getVertex(e.getTag(1))->x())/2.0, elem, k, x[0], xLoc);
        // for(int i = 0; i < _dimRecovery; ++i) monomials[i] = pow(xLoc, _expXRecovery[i]);
        for(int i = 0; i < _dimRecovery; ++i) monomials[i] = myPow(xLoc, _expXRecovery[i]);

        // std::cout<<"Monomes : "<<std::endl;
        // for(auto val : monomials)
        //   std::cout<<val<<std::endl;
        // std::cout<<std::endl;

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol) {
              case 1:
                myMat3(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
              case 2:
                myMat4(i, j) += jac * w[k] * monomials[i] * monomials[j];
                myMatContrainte(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
            }
          }
          // myMatContrainte(_dimRecovery,i) += jac * w[k] * monomials[i];
          // myMatContrainte(i,_dimRecovery) += jac * w[k] * monomials[i];
        }
      }
      // if(isBoundary(e)){
      // Not an integral :
      myMatContrainte(_dimRecovery, 0) = 1.0;
      myMatContrainte(0, _dimRecovery) = 1.0;
      // }
    }
    switch(_degSol) {
      case 1:
        lsInvAtVertices3OnEdges[e.getTag()] = myMat3.inverse();
        break;
      case 2:
        lsInvAtVertices4OnEdges[e.getTag()] = myMat4.inverse();
        lsInvAtVerticesContOnEdges[e.getTag()] = myMatContrainte.inverse();
        // std::cout<<"Edge:"<<std::endl; std::cout<<myMat4<<std::endl;
        break;
    }
  }
  printf("Least square matrices at edges : done\n");
}

void feRecovery::matrixInverseEigen2D()
{
  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord(9, 0.), x(3, 0.0), xLoc(3, 0.0), monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  std::vector<double> &solVec = _sol->getSolutionReference();

  // Matrices defined on the vertices
  printf("Looping over %ld vertices... ", _patch->getVertices().size());
  tic();
  for(auto v : _patch->getVertices()) {
    double xv = _mesh->getVertex(v)->x();
    double yv = _mesh->getVertex(v)->y();
    double zv = _mesh->getVertex(v)->z();

    Eigen::Matrix<double, 6, 6> myMat6 = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 10, 10> myMat10 = Eigen::Matrix<double, 10, 10>::Zero();
    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      // std::cout<<elem<<std::endl;
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }

      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // // Sommets de l'element
      // double r1[3] = {0., 0., 0.};
      // double r2[3] = {1., 0., 0.};
      // double r3[3] = {0., 1., 0.};
      // std::vector<double> xv1(3, 0.0);
      // std::vector<double> xv2(3, 0.0);
      // std::vector<double> xv3(3, 0.0);
      // _geoSpace->interpolateVectorField(geoCoord, r1, xv1);
      // _geoSpace->interpolateVectorField(geoCoord, r2, xv2);
      // _geoSpace->interpolateVectorField(geoCoord, r3, xv3);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        // xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        // xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        // xLoc[2] = x[2] - _mesh->getVertex(v)->z();
        xLoc[0] = (x[0] - xv);
        xLoc[1] = (x[1] - yv);
        xLoc[2] = (x[2] - zv);

        for(int i = 0; i < _dimRecovery; ++i) {
          // monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
          monomials[i] = myPow(xLoc[0], _expXRecovery[i]) * myPow(xLoc[1], _expYRecovery[i]);
        }

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol) {
              case 1:
                myMat6(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
              case 2:
                myMat10(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
            }
          }
        }
      }
    }
    switch(_degSol) {
      case 1:
        lsInvAtVertices6[v] = myMat6.inverse();
        break;
      case 2:
        lsInvAtVertices10[v] = myMat10.inverse();
        break;
    }
  }
  printf("Least square matrices at vertices : done\n");
  toc();

  // Matrices defined on the edges
  printf("Looping over %ld edges... ", _mesh->_edges.size());
  tic();
  for(auto e : _mesh->_edges) {
    // TODO : boucler sur le nombre de DOFS par edge, ici on suppose juste un P2 avec 1 dof

    Eigen::Matrix<double, 6, 6> myMat6 = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 10, 10> myMat10 = Eigen::Matrix<double, 10, 10>::Zero();

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }
      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        // printf("%d : (%f,%f) et %d : (%f,%f)\n",
        //   e.getTag(0), _mesh->getVertex(e.getTag(0))->x(), _mesh->getVertex(e.getTag(0))->y(),
        //   e.getTag(1), _mesh->getVertex(e.getTag(1))->x(), _mesh->getVertex(e.getTag(1))->y());
        xLoc[0] = x[0] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->x() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->x()) /
                           2.0;
        xLoc[1] = x[1] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->y() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->y()) /
                           2.0;
        xLoc[2] = x[2] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->z() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->z()) /
                           2.0;

        for(int i = 0; i < _dimRecovery; ++i){
          // monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
          monomials[i] = myPow(xLoc[0], _expXRecovery[i]) * myPow(xLoc[1], _expYRecovery[i]);
        }

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol) {
              case 1:
                myMat6(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
              case 2:
                myMat10(i, j) += jac * w[k] * monomials[i] * monomials[j];
                break;
            }
          }
        }
      }
    }
    switch(_degSol) {
      case 1:
        lsInvAtVertices6OnEdges[e.getTag()] = myMat6.inverse();
        break;
      case 2:
        lsInvAtVertices10OnEdges[e.getTag()] = myMat10.inverse();
        break;
    }
  }
  printf("Least square matrices at edges : done\n");
  toc();
}

void feRecovery::solveLeastSquareEigen1D(int indRecovery, int iDerivative)
{
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();
  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.0);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  std::vector<int> &vertices = _patch->getVertices();

  std::vector<double> u(indRecovery, 0.);

  std::vector<double> &solVec = _sol->getSolutionReference();

  for(auto v : vertices) {
    double xv = _mesh->getVertex(v)->x();

    std::vector<Eigen::VectorXd> RHS3(indRecovery, Eigen::MatrixXd::Zero(3, 1));
    std::vector<Eigen::VectorXd> RHS4(indRecovery, Eigen::MatrixXd::Zero(4, 1));

    std::vector<Eigen::VectorXd> RHSCont(indRecovery, Eigen::MatrixXd::Zero(5, 1));

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }

      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // Loop over quad points and increment right hand side
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        double xLoc = x[0] - xv;

        for(int i = 0; i < _dimRecovery; ++i) {
          // monomials[i] = pow(xLoc, _expXRecovery[i]);
          monomials[i] = myPow(xLoc, _expXRecovery[i]);
        }

        // printf("Vertex %d (%2.2f)- Elem %d - Quad node %d (%2.2f) - xLoc = %2.2f\n",
        //   v, xv, elem, k, x[0], xLoc);

        // std::cout<<"Monomes : "<<std::endl;
        // for(auto val : monomials)
        //   std::cout<<val<<std::endl;
        // std::cout<<std::endl;

        if(iDerivative > 0) {
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = 0.0;
          }
          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            // Vertices
            for(int iVert = 0; iVert < 2; ++iVert) {
              int vNode = _cnc->getNodeConnectivity(elem, iVert);
              std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
            }
            // Edge
            int edge = _cnc->getEdgeConnectivity(elem, 0);
            std::vector<double> &du = derivativeCoeffOnEdges[edge][0][iDeriv];
            u[iDeriv] += _intSpace->getFunctionAtQuadNode(2, k) *
                         du[0]; // FIXME : Évaluer proprement la fonction d'interpolation
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = _intSpace->interpolateFieldAtQuadNode(_solution, k);
            // u[iDeriv] = pow(x[0],5);
            // printf("u = %+-4.4f - uRef = %+-4.4f\n", u[iDeriv], pow(x[0],4));
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
          for(int i = 0; i < _dimRecovery; ++i) {
            switch(_degSol) {
              case 1:
                RHS3[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
              case 2:
                RHS4[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                RHSCont[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                // std::cout<<"RHS4 - "<<iDeriv<<std::endl;
                // std::cout<<RHS4[iDeriv]<<std::endl;
                break;
            }
          }
          // RHSCont[iDeriv](_dimRecovery) += jac * w[k] * u[iDeriv];
        }
      }
      // Not an integral:
      if(iDerivative < 3 && isBoundary(_mesh, v)) {
        if(iDerivative == 0) RHSCont[0](_dimRecovery) = boundaryCondition(_mesh, v, _solRef);
        if(iDerivative == 1) {
          std::vector<double> res(_dim, 0.0);
          boundaryConditionVec(_mesh, v, _solRefGrad, res);
          RHSCont[0](_dimRecovery) = res[0];
        }
        if(iDerivative == 2) {
          std::vector<double> res(_dim, 0.0);
          boundaryConditionVec(_mesh, v, _solRefHess, res);
          RHSCont[0](_dimRecovery) = res[0];
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
      switch(_degSol) {
        case 1: {
          Eigen::VectorXd sol = lsInvAtVertices3[v] * RHS3[iDeriv];
          recoveryCoeff[v][iDeriv] =
            std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2: {
          Eigen::VectorXd sol = lsInvAtVertices4[v] * RHS4[iDeriv];
          Eigen::VectorXd solCont = lsInvAtVerticesCont[v] * RHSCont[iDeriv];
          std::cout << "Printing system:" << std::endl;
          std::cout << lsInvAtVerticesCont[v] << std::endl;
          std::cout << std::endl;
          std::cout << RHSCont[iDeriv] << std::endl;
          std::cout << std::endl;
          std::cout << sol << std::endl;
          std::cout << std::endl;
          std::cout << solCont << std::endl;
          std::cout << std::endl;
          std::cout << std::endl;
          // recoveryCoeff[v][iDeriv] = std::vector<double>(sol.data(), sol.data() + sol.rows() *
          // sol.cols());
          std::vector<double> solContVec =
            std::vector<double>(solCont.data(), solCont.data() + solCont.rows() * solCont.cols());
          solContVec.pop_back();
          if(iDerivative < 3 && isBoundary(_mesh, v)) {
            recoveryCoeff[v][iDeriv] = solContVec;
          } else {
            recoveryCoeff[v][iDeriv] =
              std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          }
          // std::cout<<"Vertex system:"<<std::endl;
          // std::cout<<lsInvAtVertices4[v]<<std::endl;
          // std::cout<<std::endl;
          // std::cout<<RHS4[iDeriv]<<std::endl;
          // std::cout<<std::endl;
          // std::cout<<sol<<std::endl;
          // std::cout<<std::endl;
          break;
        }
      }
    }
  }

  for(auto e : _mesh->_edges) {
    // std::cout<<"Edge "<<e.getTag()<<"======================"<<std::endl;

    std::vector<Eigen::VectorXd> RHS3(indRecovery, Eigen::MatrixXd::Zero(3, 1));
    std::vector<Eigen::VectorXd> RHS4(indRecovery, Eigen::MatrixXd::Zero(4, 1));

    std::vector<Eigen::VectorXd> RHSCont(indRecovery, Eigen::MatrixXd::Zero(5, 1));

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }

      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // Loop over quad points and increment right hand side
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        double xLoc =
          x[0] - (_mesh->getVertex(e.getTag(0))->x() + _mesh->getVertex(e.getTag(1))->x()) / 2.0;

        for(int i = 0; i < _dimRecovery; ++i) {
          // monomials[i] = pow(xLoc, _expXRecovery[i]);
          monomials[i] = myPow(xLoc, _expXRecovery[i]);
        }

        if(iDerivative > 0) {
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = 0.0;
          }

          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            // Vertices
            for(int iVert = 0; iVert < 2; ++iVert) {
              int vNode = _cnc->getNodeConnectivity(elem, iVert);
              std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
            }
            // Edge
            int edge = _cnc->getEdgeConnectivity(elem, 0);
            std::vector<double> &du = derivativeCoeffOnEdges[edge][0][iDeriv];
            u[iDeriv] += _intSpace->getFunctionAtQuadNode(2, k) *
                         du[0]; // FIXME : Évaluer proprement la fonction d'interpolation
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = _intSpace->interpolateFieldAtQuadNode(_solution, k);
            // u[iDeriv] = pow(x[0],5);
            // printf("u = %+-4.4f - uRef = %+-4.4f\n", u[iDeriv], pow(x[0],4));
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
          for(int i = 0; i < _dimRecovery; ++i) {
            switch(_degSol) {
              case 1:
                RHS3[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
              case 2:
                RHS4[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                RHSCont[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                // std::cout<<"RHS4 :"<<std::endl; std::cout<<RHS4[iDeriv]<<std::endl;
                break;
            }
          }
          // RHSCont[iDeriv](_dimRecovery) += jac * w[k] * u[iDeriv];
        }
      }
      if(iDerivative < 3 && isBoundary(_mesh, e)) {
        if(iDerivative == 0) {
          RHSCont[0](_dimRecovery) = boundaryCondition(_mesh, e, _solRef);
        }
        if(iDerivative == 1) {
          std::vector<double> res(_dim, 0.0);
          boundaryConditionVec(_mesh, e, _solRefGrad, res);
          RHSCont[0](_dimRecovery) = res[0];
        }
        if(iDerivative == 2) {
          std::vector<double> res(_dim, 0.0);
          boundaryConditionVec(_mesh, e, _solRefHess, res);
          RHSCont[0](_dimRecovery) = res[0];
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
      switch(_degSol) {
        case 1: {
          Eigen::VectorXd sol = lsInvAtVertices3OnEdges[e.getTag()] * RHS3[iDeriv];
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] =
            std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2: {
          Eigen::VectorXd sol = lsInvAtVertices4OnEdges[e.getTag()] * RHS4[iDeriv];
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] =
            std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          // Eigen::VectorXd solCont = lsInvAtVerticesContOnEdges[e.getTag()]*RHSCont[iDeriv];
          // std::vector<double> solContVec = std::vector<double>(solCont.data(), solCont.data() +
          // solCont.rows() * solCont.cols()); solContVec.pop_back();
          // recoveryCoeffOnEdges[e.getTag()][0][iDeriv] = solContVec;
          // std::cout<<"Edge system:"<<std::endl;
          // std::cout<<lsInvAtVertices4OnEdges[e.getTag()]<<std::endl;
          // std::cout<<RHS4[iDeriv]<<std::endl;
          // std::cout<<"Edge sol:"<<std::endl;
          // std::cout<<sol<<std::endl;
          break;
        }
      }
    }
  }
}

void feRecovery::solveLeastSquareEigen2D(int indRecovery, int iDerivative)
{
  printf("Calling solveLS with %d\n", indRecovery);

  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.0);
  std::vector<double> xLoc(3, 0.0);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  std::vector<int> &vertices = _patch->getVertices();

  std::vector<double> u(indRecovery, 0.);

  int _nEdgePerElm = _cnc->getNbEdgePerElem();
  int _nVertPerElm = _nNodePerElm;

  std::vector<double> &solVec = _sol->getSolutionReference();

  std::vector<Eigen::VectorXd> RHS6(indRecovery, Eigen::VectorXd::Zero(6));
  std::vector<Eigen::VectorXd> RHS10(indRecovery, Eigen::VectorXd::Zero(10));

  std::vector<double> recoveryVector(10, 0.);
  Eigen::VectorXd sol;

  for(auto v : vertices) {
    double xv = _mesh->getVertex(v)->x();
    double yv = _mesh->getVertex(v)->y();
    double zv = _mesh->getVertex(v)->z();

    // std::vector<Eigen::VectorXd> RHS6(indRecovery, Eigen::MatrixXd::Zero(6, 1));
    // std::vector<Eigen::VectorXd> RHS10(indRecovery, Eigen::MatrixXd::Zero(10, 1));
    for(int i = 0; i < indRecovery; ++i){
    //   RHS6[i] = Eigen::VectorXd::Zero(6);
      RHS6[i].setZero(6);
    //   RHS10[i] = Eigen::VectorXd::Zero(10);
      RHS10[i].setZero(10);
    }

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }
      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // // Sommets de l'element
      // double r1[3] = {0., 0., 0.};
      // double r2[3] = {1., 0., 0.};
      // double r3[3] = {0., 1., 0.};
      // std::vector<double> xv1(3, 0.0);
      // std::vector<double> xv2(3, 0.0);
      // std::vector<double> xv3(3, 0.0);
      // _geoSpace->interpolateVectorField(geoCoord, r1, xv1);
      // _geoSpace->interpolateVectorField(geoCoord, r2, xv2);
      // _geoSpace->interpolateVectorField(geoCoord, r3, xv3);

      // Loop over quad points and increment right hand side
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        // xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        // xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        // xLoc[2] = x[2] - _mesh->getVertex(v)->z();
        xLoc[0] = (x[0] - xv);
        xLoc[1] = (x[1] - yv);
        xLoc[2] = (x[2] - zv);

        for(int i = 0; i < _dimRecovery; ++i) {
          // monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
          monomials[i] = myPow(xLoc[0], _expXRecovery[i]) * myPow(xLoc[1], _expYRecovery[i]);
        }

        if(iDerivative > 0) {
          // The contributions to the derivative from the vertices must be evaluated and averaged at
          // quad nodes to avoid a trivial solution
          // u = 0.0;
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = 0.0;
          }
          // for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
          //   std::cout<<_nNodePerElm<<std::endl;
          //   int vNode = _cnc->getNodeConnectivity(elem, iNode);
          //   // Get the coefficients of the derivative (used only if recovering a derivative)
          //   for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
          //     std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
          //     // std::vector<double> &du = derivativeCoeff[vNode][indRecovery];
          //     xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
          //     xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
          //     xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
          //     for(int i = 0; i < _dimDerivation; ++i) {
          //       u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0],
          //       _expX[i]) *
          //            pow(xLoc[1], _expY[i]);
          //     }
          //   }
          // }

          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            // Vertices
            // printf("elem with vertices\n");
            for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
              int vNode = _cnc->getNodeConnectivity(elem, iVert);
              std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
              // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
              // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
              // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
              // printf("%d en %f - %f\n", vNode, _mesh->getVertex(vNode)->x(),
              // _mesh->getVertex(vNode)->y()); for(int i = 0; i < _dimDerivation; ++i) { u[iDeriv]
              // += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) *
              // pow(xLoc[1], _expY[i]); u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) *
              // du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              // }
            }
            // Edges
            for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
              int localEdge = fabs(_cnc->getEdgeConnectivity(elem, iEdge));
              std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
              // int v1, v2;
              // if(iEdge == _nEdgePerElm-1){
              //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
              //   v2 = _cnc->getNodeConnectivity(elem, 0);
              // } else{
              //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
              //   v2 = _cnc->getNodeConnectivity(elem, iEdge+1);
              // }
              // xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
              // xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
              // xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
              // printf("edge %d - %d en %f - %f et %f - %f\n",
              //   v1, v2, _mesh->getVertex(v1)->x(), _mesh->getVertex(v1)->y(),
              //   _mesh->getVertex(v2)->x(), _mesh->getVertex(v2)->y());

              // Edge e(_mesh->getVertex(v1),_mesh->getVertex(v2));
              // std::set<Edge, EdgeLessThan>::iterator ret;
              // ret = _mesh->_edges.find(e);
              // if(ret != _mesh->_edges.end()){
              //   printf("sanityCheck : edge %d was found at %d : relie les sommets %d et %d\n",
              //   localEdge, ret->getTag(), ret->getTag(0), ret->getTag(1));
              // } else{
              //   printf("sanityCheck : edge %d was NOT found\n", localEdge, _mesh->_edges);
              // }
              // for(int i = 0; i < _dimDerivation; ++i) {
              // Attention : iEdge+3 hardcoded for P2 interpolant
              // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
              // _expX[i]) * pow(xLoc[1], _expY[i]); u[iDeriv] +=
              // _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) *
              // pow(xLoc[1], _expY[i]);
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) * du[0];
              // }
            }
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = _intSpace->interpolateFieldAtQuadNode(_solution, k);
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
          // u[iDeriv] = 2.0;
          for(int i = 0; i < _dimRecovery; ++i) {
            switch(_degSol) {
              case 1:
                RHS6[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
              case 2:
                RHS10[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
            }
          }
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
      switch(_degSol) {
        case 1: {
          sol = lsInvAtVertices6[v] * RHS6[iDeriv];
          recoveryCoeff[v][iDeriv] = std::vector<double>(sol.data(), sol.data() + 6);
            // std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2: {
          sol = lsInvAtVertices10[v] * RHS10[iDeriv];
          for(int i = 0; i < 10; ++i){
            recoveryVector[i] = sol(i);
          }
          recoveryCoeff[v][iDeriv] = recoveryVector;
          // recoveryCoeff[v][iDeriv] =
          //   std::vector<double>(sol.data(), sol.data() + 10);
            // std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
      }
    }
  }

  for(auto e : _mesh->_edges) {
    // std::vector<Eigen::VectorXd> RHS6(indRecovery, Eigen::MatrixXd::Zero(6, 1));
    // std::vector<Eigen::VectorXd> RHS10(indRecovery, Eigen::MatrixXd::Zero(10, 1));
    for(int i = 0; i < indRecovery; ++i){
    //   RHS6[i] = Eigen::VectorXd::Zero(6);
      RHS6[i].setZero(6);
    //   RHS10[i] = Eigen::VectorXd::Zero(10);
      RHS10[i].setZero(10);
    }

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem, _adr);
      for(size_t i = 0; i < _adr.size(); ++i) {
        _solution[i] = solVec[_adr[i]];
      }
      _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

      // // Sommets de l'element
      // double r1[3] = {0., 0., 0.};
      // double r2[3] = {1., 0., 0.};
      // double r3[3] = {0., 1., 0.};
      // std::vector<double> xv1(3, 0.0);
      // std::vector<double> xv2(3, 0.0);
      // std::vector<double> xv3(3, 0.0);
      // _geoSpace->interpolateVectorField(geoCoord, r1, xv1);
      // _geoSpace->interpolateVectorField(geoCoord, r2, xv2);
      // _geoSpace->interpolateVectorField(geoCoord, r3, xv3);

      // Loop over quad points and increment right hand side
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        xLoc[0] = x[0] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->x() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->x()) /
                           2.0;
        xLoc[1] = x[1] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->y() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->y()) /
                           2.0;
        xLoc[2] = x[2] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->z() +
                          _mesh->getVertexFromGmshNodeTag(e.getTag(1))->z()) /
                           2.0;

        for(int i = 0; i < _dimRecovery; ++i) {
          // monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
          monomials[i] = myPow(xLoc[0], _expXRecovery[i]) * myPow(xLoc[1], _expYRecovery[i]);
        }

        if(iDerivative > 0) {
          // The contributions to the derivative from the vertices must be evaluated and averaged at
          // quad nodes to avoid a trivial solution
          // u = 0.0;
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = 0.0;
          }
          // for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
          //   int vNode = _cnc->getNodeConnectivity(elem, iNode);
          //   // Get the coefficients of the derivative (used only if recovering a derivative)
          //   for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
          //     std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
          //     // std::vector<double> &du = derivativeCoeff[vNode][indRecovery];
          //     xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
          //     xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
          //     xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
          //     for(int i = 0; i < _dimDerivation; ++i) {
          //       u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0],
          //       _expX[i]) * pow(xLoc[1], _expY[i]);
          //     }
          //   }
          // }

          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            // Vertices
            // printf("elem with vertices\n");
            for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
              int vNode = _cnc->getNodeConnectivity(elem, iVert);

              // OPT_CHANGE
              // std::vector<double> &du = derivativeCoeff[vNode][iDeriv];

              // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
              // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
              // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
              // printf("%d en %f - %f\n", vNode, _mesh->getVertex(vNode)->x(),
              // _mesh->getVertex(vNode)->y()); for(int i = 0; i < _dimDerivation; ++i) { u[iDeriv]
              // += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) *
              // pow(xLoc[1], _expY[i]); u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) *
              // du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);

              // OPT_CHANGE
              // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][iDeriv][0];

              // }
            }
            // Edges
            for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
              int localEdge = fabs(_cnc->getEdgeConnectivity(elem, iEdge));

              // OPT_CHANGE
              // std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
              // du = derivativeCoeffOnEdges[localEdge][0][iDeriv];

              // int v1, v2;
              // if(iEdge == _nEdgePerElm-1){
              //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
              //   v2 = _cnc->getNodeConnectivity(elem, 0);
              // } else{
              //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
              //   v2 = _cnc->getNodeConnectivity(elem, iEdge+1);
              // }
              // xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
              // xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
              // xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
              // printf("edge %d - %d en %f - %f et %f - %f\n",
              //   v1, v2, _mesh->getVertex(v1)->x(), _mesh->getVertex(v1)->y(),
              //   _mesh->getVertex(v2)->x(), _mesh->getVertex(v2)->y());

              // Edge e(_mesh->getVertex(v1),_mesh->getVertex(v2));
              // std::set<Edge, EdgeLessThan>::iterator ret;
              // ret = _mesh->_edges.find(e);
              // if(ret != _mesh->_edges.end()){
              //   printf("sanityCheck : edge %d was found at %d : relie les sommets %d et %d\n",
              //   localEdge, ret->getTag(), ret->getTag(0), ret->getTag(1));
              // } else{
              //   printf("sanityCheck : edge %d was NOT found\n", localEdge, _mesh->_edges);
              // }
              // for(int i = 0; i < _dimDerivation; ++i) {
              // Attention : iEdge+3 hardcoded for P2 interpolant
              // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
              // _expX[i]) * pow(xLoc[1], _expY[i]); u[iDeriv] +=
              // _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) *
              // pow(xLoc[1], _expY[i]);

              // OPT_CHANGE
              // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) * du[0];
              u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) * derivativeCoeffOnEdges[localEdge][0][iDeriv][0];
              // }
            }
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
            u[iDeriv] = _intSpace->interpolateFieldAtQuadNode(_solution, k);
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
          // u[iDeriv] = 2.0;
          for(int i = 0; i < _dimRecovery; ++i) {
            switch(_degSol) {
              case 1:
                RHS6[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
              case 2:
                RHS10[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i];
                break;
            }
          }
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv) {
      switch(_degSol) {
        case 1: {
          sol = lsInvAtVertices6OnEdges[e.getTag()] * RHS6[iDeriv];
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] =
            std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2: {
          sol = lsInvAtVertices10OnEdges[e.getTag()] * RHS10[iDeriv];
          for(int i = 0; i < 10; ++i){
            recoveryVector[i] = sol(i);
          }
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] = recoveryVector;
          // recoveryCoeffOnEdges[e.getTag()][0][iDeriv] =
          //   std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
      }
    }
  }
}

// void feRecovery::solveLeastSquare(int indRecovery, bool recoverDerivative) {
//   PetscErrorCode ierr;

//   if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
//     printf("Error : mismatch in number of quadrature points\n");
//   }
//   std::vector<double> &w = _geoSpace->getQuadratureWeights();
//   std::vector<double> &J = _cnc->getJacobians();

//   std::vector<double> geoCoord;
//   std::vector<double> x(3, 0.0);
//   std::vector<double> xLoc(3, 0.0);
//   std::vector<double> monomials(_dimRecovery, 0.);

//   // La solution interpolée aux points d'intégration
//   double u, jac;

//   int nQuad = _geoSpace->getNbQuadPoints();

//   std::vector<int> &vertices = _patch->getVertices();
//   for(auto v : vertices) {
//     // std::cout<<"Assembling on vertex "<<v<<std::endl;
//     // Reset matrix and RHS
//     ierr = MatZeroEntries(A);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     ierr = VecZeroEntries(b);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     // Get patch of elements
//     std::set<int> &elemPatch = _patch->getPatch(v);

//     // std::string name = "patchesVert" + std::to_string(v) + ".pos";

//     // FILE *file = fopen(name.c_str(),"w");
//     // std::string text = "View\"patchVertex" + std::to_string(v) + "\"{\n";
//     // fprintf(file, text.c_str());

//     for(auto elem : elemPatch) {
//       // std::cout<<elem<<std::endl;
//       _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
//                                             elem,  _adr);
//       _intSpace->initializeSolution(_sol, _adr);
//       _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

//       // Sommets de l'element
//       double r1[3] = {0., 0., 0.};
//       double r2[3] = {1., 0., 0.};
//       double r3[3] = {0., 1., 0.};
//       std::vector<double> xv1(3, 0.0);
//       std::vector<double> xv2(3, 0.0);
//       std::vector<double> xv3(3, 0.0);
//       _geoSpace->interpolateVectorField(geoCoord, r1, xv1);
//       _geoSpace->interpolateVectorField(geoCoord, r2, xv2);
//       _geoSpace->interpolateVectorField(geoCoord, r3, xv3);
//       // fprintf(file,
//       //    "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n",
//       //    xv1[0], xv1[1], 0.0, xv2[0], xv2[1], 0.0, xv3[0], xv3[1], 0.0, elem, elem, elem);

//       // Loop over quad points and increment least square matrix
//       for(int k = 0; k < nQuad; ++k) {
//         jac = J[nQuad * elem + k];

//         // TODO : normaliser ?
//         _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
//         xLoc[0] = x[0] - _mesh->getVertex(v)->x();
//         xLoc[1] = x[1] - _mesh->getVertex(v)->y();
//         xLoc[2] = x[2] - _mesh->getVertex(v)->z();

//         for(int i = 0; i < _dimRecovery; ++i)
//           monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);

//         for(int i = 0; i < _dimRecovery; ++i) {
//           // valb[i] = jac * w[k] * u * monomials[i];
//           for(int j = 0; j < _dimRecovery; ++j) {
//             valA[_dimRecovery * i + j] = jac * w[k] * monomials[i] * monomials[j];
//           }
//         }
//         ierr = MatSetValues(A, _dimRecovery, indI, _dimRecovery, indJ, valA, ADD_VALUES);

//         if(recoverDerivative) {
//           // The contributions to the derivative from the vertices must be evaluated and averaged
//           at
//           // quad nodes to avoid a trivial solution
//           u = 0.0;
//           for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
//             int vNode = _cnc->getNodeConnectivity(elem, iNode);
//             // Get the coefficients of the derivative (used only if recovering a derivative)
//             std::vector<double> &du = derivativeCoeff[vNode][indRecovery];
//             // printf("iNode = %d \t vertex = %d\n", iNode, vNode);
//             xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
//             xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
//             xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
//             // printf("iVert = %d - k = %2d - iNode = %d - xInt = %+-4.4f yInt = %+-4.4f xLoc =
//             // %+-4.4f yLoc = %+-4.4f\n", v, k, iNode, x[0], x[1], xLoc[0], xLoc[1]);
//             for(int i = 0; i < _dimDerivation; ++i) {
//               u += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0], _expX[i]) *
//                    pow(xLoc[1], _expY[i]);
//               // printf("Rec %d : %+-4.4f \t %+-4.4f \t %+-4.4f \t %+-4.4f\n", indRecovery, u,
//               // _geoSpace->getFunctionAtQuadNode(iNode,k), du[i], pow(xLoc[0],
//               // _expX[i])*pow(xLoc[1], _expY[i]));
//             }
//           }
//         } else {
//           // Simply interpolate the solution at quad nodes
//           u = _intSpace->interpolateSolutionAtQuadNode(k);
//         }

//         for(int i = 0; i < _dimRecovery; ++i) valb[i] = jac * w[k] * u * monomials[i];

//         ierr = VecSetValues(b, _dimRecovery, indI, valb, ADD_VALUES);
//       }
//     }

//     // MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//     // MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//     // MatView(A, PETSC_VIEWER_STDOUT_WORLD);

//     // fprintf(file, "};\n");
//     // fclose(file);

//     ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);

//     // Solve system
//     PetscInt its;
//     ierr = KSPSolve(ksp, b, c);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//     KSPGetIterationNumber(ksp, &its);
//     // PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
//     VecSet(res, 0.0);
//     MatMult(A, c, res);
//     VecAXPY(res, -1.0, b);
//     double normeAxmb = 0.0;
//     ierr = VecNorm(res, NORM_MAX, &normeAxmb);
//     CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     // std::cout<<"Norme du résidu matriciel Ax-b : "<<normeAxmb<<std::endl;

//     // printf("Solution at vertex %f - %f - %f : \n", _mesh->getVertex(v)->x(),
//     // _mesh->getVertex(v)->y(), _mesh->getVertex(v)->z()); VecView(c,PETSC_VIEWER_STDOUT_WORLD);

//     // Copy coefficients in main structure
//     PetscScalar *array;
//     VecGetArray(c, &array);
//     std::vector<double> coeffs((double *)array, (double *)array + _dimRecovery);
//     VecRestoreArray(c, &array);
//     recoveryCoeff[v][indRecovery] = coeffs;

//     if(v == 0){
//       for(auto val : coeffs){
//         std::cout<<"Solution PETSc : "<<std::endl;
//         std::cout<<val<<std::endl;
//       }
//     }

//   }
// }

void feRecovery::derivative(int indRecovery, int iDerivative, std::ostream &output)
{
  std::vector<int> &vertices = _patch->getVertices();
  double tol = 0.0;

  if(iDerivative == 0) {
    /* Store the independent term of the solution.
    Used in feRecovery::evalDerivative. */
    std::vector<double> solV(vertices.size(), 0.0);
    for(auto v : vertices) solV[v] = recoveryCoeff[v][0][0];
    derivAtVertices.push_back(solV);
    std::vector<double> solE(_mesh->getNbEdges(), 0.0);
    for(auto e : _mesh->_edges) solE[e.getTag() - 1] = recoveryCoeffOnEdges[e.getTag()][0][0][0];
    derivAtEdges.push_back(solE);
  }

  if(_dim == 1) {
    for(auto v : vertices) {
      std::vector<double> &u = recoveryCoeff[v][indRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      for(int i = 0; i < _dimDerivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeff[v][indRecovery] = dudx;
    }

    for(auto e : _mesh->_edges) {
      std::vector<double> &u = recoveryCoeffOnEdges[e.getTag()][0][indRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      for(int i = 0; i < _dimDerivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeffOnEdges[e.getTag()][0][indRecovery] = dudx;
    }

  } else if(_dim == 2) {
    // To store the derivatives
    std::vector<double> dudxV(vertices.size(), 0.0);
    std::vector<double> dudyV(vertices.size(), 0.0);
    std::vector<double> dudxE(_mesh->getNbEdges(), 0.0);
    std::vector<double> dudyE(_mesh->getNbEdges(), 0.0);

    for(auto v : vertices) {
      std::vector<double> &u = recoveryCoeff[v][indRecovery];
      int indX = 0, indY = 0;
      std::vector<double> dudx(_dimDerivation, 0.);
      std::vector<double> dudy(_dimDerivation, 0.);
      double sumCoeffX = 0.0;
      double sumCoeffY = 0.0;
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) {
          sumCoeffX += ((double)_expXRecovery[i]) * fabs(u[i]);
        }
        if(_expYRecovery[i] != 0) {
          sumCoeffY += ((double)_expYRecovery[i]) * fabs(u[i]);
        }
      }
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) {
          if(fabs(u[i] / sumCoeffX) > tol) {
            dudx[indX] = _expXRecovery[i] * u[i];
          }
          indX++;
        }
        if(_expYRecovery[i] != 0) {
          if(fabs(u[i] / sumCoeffY) > tol) {
            dudy[indY] = _expYRecovery[i] * u[i];
          }
          indY++;
        }
      }
      derivativeCoeff[v][2 * indRecovery] = dudx;
      derivativeCoeff[v][2 * indRecovery + 1] = dudy;

      dudxV[v] = dudx[0];
      dudyV[v] = dudy[0];

      // printf("Derivee à l'emplacement %d et %d = %+-4.4f - %+-4.4f\n", 2*indRecovery,
      // 2*indRecovery+1, derivativeCoeff[v][2 * indRecovery][0], derivativeCoeff[v][2 *
      // indRecovery+1][0]);
    }

    for(auto e : _mesh->_edges) {
      int indX = 0, indY = 0;
      std::vector<double> &u = recoveryCoeffOnEdges[e.getTag()][0][indRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      std::vector<double> dudy(_dimDerivation, 0.);
      double sumCoeffX = 0.0;
      double sumCoeffY = 0.0;
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) {
          sumCoeffX += ((double)_expXRecovery[i]) * fabs(u[i]);
        }
        if(_expYRecovery[i] != 0) {
          sumCoeffY += ((double)_expYRecovery[i]) * fabs(u[i]);
        }
      }
      // for(int i = 0; i < _dimRecovery; ++i) {
      //   if(_expXRecovery[i] != 0 && fabs(u[i]/sumCoeffX) > tol) { dudx[indX++] = _expXRecovery[i]
      //   * u[i]; } if(_expYRecovery[i] != 0 && fabs(u[i]/sumCoeffY) > tol) { dudy[indY++] =
      //   _expYRecovery[i] * u[i]; }
      // }
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) {
          if(fabs(u[i] / sumCoeffX) > tol) {
            dudx[indX] = _expXRecovery[i] * u[i];
          }
          indX++;
        }
        if(_expYRecovery[i] != 0) {
          if(fabs(u[i] / sumCoeffY) > tol) {
            dudy[indY] = _expYRecovery[i] * u[i];
          }
          indY++;
        }
      }
      derivativeCoeffOnEdges[e.getTag()][0][2 * indRecovery] = dudx;
      derivativeCoeffOnEdges[e.getTag()][0][2 * indRecovery + 1] = dudy;

      dudxE[e.getTag() - 1] = dudx[0];
      dudyE[e.getTag() - 1] = dudy[0];
    }

    derivAtVertices.push_back(dudxV);
    derivAtVertices.push_back(dudyV);
    derivAtEdges.push_back(dudxE);
    derivAtEdges.push_back(dudyE);

    // Exporter les derivees aux sommets (terme indépendant)
    std::string filename = "derivees_var" + _intSpace->getFieldID() + "_ordre" +
                           std::to_string(iDerivative) + "_indRec" + std::to_string(indRecovery) +
                           ".txt";
    FILE *myfile = fopen(filename.c_str(), "w");
    printf("Writing derivatives to %s\n", filename.c_str());
    for(auto v : vertices) {
      fprintf(myfile, "%+-12.12f \t %+-12.12f\n", derivativeCoeff[v][2 * indRecovery][0],
              derivativeCoeff[v][2 * indRecovery + 1][0]);
    }
    fclose(myfile);
  } else {
    printf("TODO : Derivatives for recoveries in 0D or 3D\n");
  }

  // Print the solution at first pass
  if(iDerivative == 0) {
    output << "$NodeData\n";
    output << "1\n\"" << _intSpace->getFieldID() << "\"\n1\n3000\n3\n0\n1\n"
           << vertices.size() << "\n";
    for(auto v : vertices) {
      // printf("Writing %d at vertex %d\n", v, _mesh->getVertex(cnt)->getTag());
      // output << _mesh->getVertex(cnt++)->getTag() << " " << recoveryCoeff[v][indRecovery][0] <<
      // std::endl;
      output << _mesh->getVertex(v)->getTag() << " " << recoveryCoeff[v][indRecovery][0]
             << std::endl;
      // output << v << " " << recoveryCoeff[v][indRecovery][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }
  // Also write the derivatives to the common file
  for(int i = 0; i < _dim; ++i) {
    output << "$NodeData\n";
    std::string fieldName = "d" + std::to_string(iDerivative + 1) + _intSpace->getFieldID() +
                            suffix[{iDerivative + 1, _dim * indRecovery + i}];
    output << "1\n\"" << fieldName << "\"\n1\n3000\n3\n0\n1\n" << vertices.size() << "\n";
    for(auto v : vertices) {
      // output << _mesh->getVertex(cnt++)->getTag() << " " << derivativeCoeff[v][_dim * indRecovery
      // + i][0] << std::endl;
      output << _mesh->getVertex(v)->getTag() << " "
             << derivativeCoeff[v][_dim * indRecovery + i][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }
}

void feRecovery::secondDerivative(int indRecovery, int iDerivative, std::ostream &output)
{
  std::vector<int> &vertices = _patch->getVertices();

  if(_dim == 1) {
    for(auto v : vertices) {
      std::vector<double> &u = derivativeCoeff[v][indRecovery];
      std::vector<double> dudx(_dim2Derivation, 0.);
      for(int i = 0; i < _dim2Derivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeff[v][indRecovery] = dudx;
    }

    for(auto e : _mesh->_edges) {
      std::vector<double> &u = derivativeCoeffOnEdges[e.getTag()][0][indRecovery];
      std::vector<double> dudx(_dim2Derivation, 0.);
      for(int i = 0; i < _dim2Derivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeffOnEdges[e.getTag()][0][indRecovery] = dudx;
    }
  }
}

void feRecovery::getErrorPolynomials()
{
  if(_dim == 2) {
    std::vector<int> &vertices = _patch->getVertices();
    std::vector<double> error(_degSol + 2, 0.);
    switch(_degSol) {
      case 1:
        for(auto v : vertices) {
          error[0] = derivativeCoeff[v][0][0];
          error[1] = derivativeCoeff[v][1][0] + derivativeCoeff[v][2][0];
          error[2] = derivativeCoeff[v][3][0];
          errorCoeff[v] = error;
        }
        break;
      case 2:
        for(auto v : vertices) {
          error[0] = derivativeCoeff[v][0][0];
          error[1] = derivativeCoeff[v][1][0] + derivativeCoeff[v][2][0] + derivativeCoeff[v][4][0];
          error[2] = derivativeCoeff[v][3][0] + derivativeCoeff[v][5][0] + derivativeCoeff[v][6][0];
          error[3] = derivativeCoeff[v][7][0];
          errorCoeff[v] = error;
        }
        break;
      case 3:
        for(auto v : vertices) {
          error[0] = derivativeCoeff[v][0][0];
          error[1] = derivativeCoeff[v][1][0] + derivativeCoeff[v][2][0] +
                     derivativeCoeff[v][4][0] + derivativeCoeff[v][8][0];
          error[2] = derivativeCoeff[v][3][0] + derivativeCoeff[v][5][0] +
                     derivativeCoeff[v][6][0] + derivativeCoeff[v][9][0] +
                     derivativeCoeff[v][10][0] + derivativeCoeff[v][12][0];
          error[3] = derivativeCoeff[v][7][0] + derivativeCoeff[v][11][0] +
                     derivativeCoeff[v][13][0] + derivativeCoeff[v][14][0];
          error[4] = derivativeCoeff[v][15][0];
          errorCoeff[v] = error;
        }
        break;
      default:
        printf("Error : Computation of error coefficients is not implemented for this solution "
               "polynomial degree (%d).\n",
               _degSol);
    }
  } else {
    printf("TODO : Compute error coefficients in 1D and 3D\n");
  }
}

void feRecovery::estimateError(std::vector<double> &norm, feFunction *solRef)
{
  norm[0] = norm[1] = 0.0;

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  std::vector<double> &solVec = _sol->getSolutionReference();

  FILE *f = fopen("solutionReconstruite.txt", "w");

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    // for(auto iElm : elemPatch) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }

    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // u reconstruit au point d'intégration : interpolation des valeurs issues des sommets
      double uReconstruit = 0.0;
      // for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
      //   int v = _cnc->getNodeConnectivity(iElm, iNode);

      //   xLoc[0] = x[0] - _mesh->getVertex(v)->x();
      //   xLoc[1] = x[1] - _mesh->getVertex(v)->y();
      //   xLoc[2] = x[2] - _mesh->getVertex(v)->z();

      //   for(int i = 0; i < _dimRecovery; ++i) {
      //     monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
      //     uReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * recoveryCoeff[v][0][i] *
      //     monomials[i];
      //   }
      // }

      // Good way ?
      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // printf("vert %d of %d\n", iVert, _nVertPerElm);
        // printf("u = %4.4f\n", recoveryCoeff[vNode][0][0]);
        uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][0];
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // int localEdge = _cnc->getEdgeConnectivity(iElm, 0);
        // printf("edge %d - %d of %d\n", iEdge, localEdge, _nEdgePerElm);
        // printf("u = %4.4f\n", recoveryCoeffOnEdges[localEdge][0][0][0]);
        if(_dim == 1)
          uReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 2, k) *
                          recoveryCoeffOnEdges[localEdge][0][0][0];
        if(_dim == 2)
          uReconstruit +=
            _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
            recoveryCoeffOnEdges[localEdge][0][0][0]; // FIXME : make this clean (iEdge+2)
      }

      // // Test
      // // Get the coefficients of the derivative (used only if recovering a derivative)
      // // Vertices
      // for(int iVert = 0; iVert < _nVertPerElm; ++iVert){
      //   int vNode = _cnc->getNodeConnectivity(iElm, iVert);
      //   xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
      //   // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
      //   // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
      //   // printf("%d en %f - %f\n", vNode, _mesh->getVertex(vNode)->x(),
      //   _mesh->getVertex(vNode)->y()); for(int i = 0; i < _dimDerivation; ++i) {
      //     // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0],
      //     _expX[i]) * pow(xLoc[1], _expY[i]);
      //     // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0],
      //     _expX[i]) * pow(xLoc[1], _expY[i]);
      //     // uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
      //     recoveryCoeff[vNode][0][i] * pow(xLoc[0], _expX[i]); uReconstruit +=
      //     _geoSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][i] * pow(xLoc[0],
      //     _expX[i]);
      //   }
      // }

      // // Edges
      // for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
      //   int localEdge = fabs(_cnc->getEdgeConnectivity(elem, iEdge));
      //   std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
      //   // int v1, v2;
      //   // if(iEdge == _nEdgePerElm-1){
      //   //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
      //   //   v2 = _cnc->getNodeConnectivity(elem, 0);
      //   // } else{
      //   //   v1 = _cnc->getNodeConnectivity(elem, iEdge);
      //   //   v2 = _cnc->getNodeConnectivity(elem, iEdge+1);
      //   // }
      //   // xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
      //   // xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
      //   // xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
      //   // printf("edge %d - %d en %f - %f et %f - %f\n",
      //   //   v1, v2, _mesh->getVertex(v1)->x(), _mesh->getVertex(v1)->y(),
      //   _mesh->getVertex(v2)->x(), _mesh->getVertex(v2)->y());

      //   // Edge e(_mesh->getVertex(v1),_mesh->getVertex(v2));
      //   // std::set<Edge, EdgeLessThan>::iterator ret;
      //   // ret = _mesh->_edges.find(e);
      //   // if(ret != _mesh->_edges.end()){
      //   //   printf("sanityCheck : edge %d was found at %d : relie les sommets %d et %d\n",
      //   localEdge, ret->getTag(), ret->getTag(0), ret->getTag(1));
      //   // } else{
      //   //   printf("sanityCheck : edge %d was NOT found\n", localEdge, _mesh->_edges);
      //   // }
      //   // for(int i = 0; i < _dimDerivation; ++i) {
      //     // Attention : iEdge+3 hardcoded for P2 interpolant
      //     // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
      //     _expX[i]) * pow(xLoc[1], _expY[i]);
      //     // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
      //     _expX[i]) * pow(xLoc[1], _expY[i]); u[iDeriv] +=
      //     _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[0];
      //   // }
      // }

      norm[0] += J[nQuad * iElm + k] * w[k] *
                 pow(uReconstruit - _intSpace->interpolateFieldAtQuadNode(_solution, k), 2);
      if(solRef) norm[1] += J[nQuad * iElm + k] * w[k] * pow(uReconstruit - solRef->eval(0, x), 2);

      fprintf(f, "%+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e\n", x[0], x[1],
              uReconstruit, solRef->eval(0, x),
              _intSpace->interpolateFieldAtQuadNode(_solution, k));
    }
  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);

  fclose(f);
}

void feRecovery::estimateDudxError(std::vector<double> &norm, feVectorFunction *solRefGrad)
{
  norm[14] = norm[15] = norm[16] = norm[17] = 0.0;

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  // for(int iElm = 0; iElm < _nElm; ++iElm) {

  std::set<int> &elemPatch = _patch->getPatch(4);

  std::vector<double> &solVec = _sol->getSolutionReference();

  // for(int iElm = 0; iElm < _nElm; ++iElm) {
  for(auto iElm : elemPatch) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }

    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // u reconstruit au point d'intégration : interpolation des valeurs issues des sommets
      double dudxReconstruit = 0.0;
      double dudyReconstruit = 0.0;

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
        // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][i] *
        // monomials[i];
        dudxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][0];
        if(_dim == 2)
          dudyReconstruit +=
            _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][1][0];
        // }
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // int v1, v2;
        // if(iEdge == _nEdgePerElm-1){
        //   v1 = _cnc->getNodeConnectivity(iElm, iEdge);
        //   v2 = _cnc->getNodeConnectivity(iElm, 0);
        // } else{
        //   v1 = _cnc->getNodeConnectivity(iElm, iEdge);
        //   v2 = _cnc->getNodeConnectivity(iElm, iEdge+1);
        // }
        // xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
        // xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
        // xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
        // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // uReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        // recoveryCoeffOnEdges[localEdge][0][0][i] * monomials[i];
        if(_dim == 1)
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 2, k) *
                             recoveryCoeffOnEdges[localEdge][0][0][0];

        if(_dim == 2)
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                             recoveryCoeffOnEdges[localEdge][0][0][0];
        if(_dim == 2)
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                             recoveryCoeffOnEdges[localEdge][0][1][0];
        // }
      }

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      double jac = 0.0, duhdx = 0.0, duhdy = 0.0;

      switch(_dim) {
        case 1:
          jac = dxdr[0];
          duhdx = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k);
          duhdx /= jac;
          break;
        case 2:
          jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
          double drdx = dxds[1] / jac;
          double drdy = -dxds[0] / jac;
          double dsdx = -dxdr[1] / jac;
          double dsdy = dxdr[0] / jac;
          duhdx = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k) * drdx +
                  _intSpace->interpolateFieldAtQuadNode_sDerivative(_solution, k) * dsdx;
          duhdy = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k) * drdy +
                  _intSpace->interpolateFieldAtQuadNode_sDerivative(_solution, k) * dsdy;
          break;
      }

      norm[14] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - duhdx, 2);
      norm[15] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - duhdy, 2);
      if(solRefGrad) {
        std::vector<double> res(2, 0.);
        solRefGrad->eval(0, x, res);
        norm[16] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - res[0], 2);
        norm[17] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - res[1], 2);
      }
    }
  }
  norm[14] = sqrt(norm[14]);
  norm[15] = sqrt(norm[15]);
  norm[16] = sqrt(norm[16]);
  norm[17] = sqrt(norm[17]);
}

void feRecovery::estimateH1Error(std::vector<double> &norm, feVectorFunction *solRefGrad)
{
  norm[2] = norm[3] = norm[4] = norm[5] = norm[6] = norm[7] = 0.0;
  norm[18] = norm[19] = norm[20] = 0.0;

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimDerivation, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  FILE *f = fopen("duReconstruite.txt", "w");

  // for(int iElm = 0; iElm < _nElm; ++iElm) {

  std::vector<double> &solVec = _sol->getSolutionReference();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    // for(auto iElm : elemPatch) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }

    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // derivee reconstruite au point d'intégration : interpolation des valeurs issues des sommets
      double dudxReconstruit = 0.0;
      double dudyReconstruit = 0.0;
      // for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
      //   int v = _cnc->getNodeConnectivity(iElm, iNode);

      //   xLoc[0] = x[0] - _mesh->getVertex(v)->x();
      //   xLoc[1] = x[1] - _mesh->getVertex(v)->y();
      //   xLoc[2] = x[2] - _mesh->getVertex(v)->z();

      //   for(int i = 0; i < _dimDerivation; ++i) {
      //     monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
      //     dudxReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) *
      //     derivativeCoeff[v][0][i] * monomials[i]; dudyReconstruit +=
      //     _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][1][i] * monomials[i];
      //   }
      // }

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // std::vector<double> &dux = derivativeCoeff[vNode][0];
        // std::vector<double> &duy = derivativeCoeff[vNode][1];
        xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
        // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // dudxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
        // derivativeCoeff[vNode][0][i] * monomials[i];
        dudxReconstruit +=
          _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][0];
        // dudyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
        // derivativeCoeff[vNode][1][i] * monomials[i];
        if(_dim == 2)
          dudyReconstruit +=
            _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][0];
        // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i])
        // * pow(xLoc[1], _expY[i]);
        // }
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
        int v1, v2;
        if(iEdge == _nEdgePerElm - 1) {
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, 0);
        } else {
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, iEdge + 1);
        }
        xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x()) / 2.0;
        xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y()) / 2.0;
        xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z()) / 2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
        // Attention : iEdge+3 hardcoded for P2 interpolant
        // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
        // _expX[i]) * pow(xLoc[1], _expY[i]); u[iDeriv] +=
        // _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) *
        // pow(xLoc[1], _expY[i]); monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        // derivativeCoeffOnEdges[localEdge][0][0][i] * monomials[i];
        if(_dim == 1)
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 2, k) *
                             derivativeCoeffOnEdges[localEdge][0][0][0];

        if(_dim == 2)
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                             derivativeCoeffOnEdges[localEdge][0][0][0];
        // dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        // derivativeCoeffOnEdges[localEdge][0][1][i] * monomials[i];
        if(_dim == 2)
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                             derivativeCoeffOnEdges[localEdge][0][1][0];
        // }
      }

      // std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      // std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      // _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      // _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      // double jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      // double drdx = dxds[1] / jac;
      // double drdy = -dxds[0] / jac;
      // double dsdx = -dxdr[1] / jac;
      // double dsdy = dxdr[0] / jac;

      // double duhdx = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
      //                _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      // double duhdy = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
      //                _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      double jac = 0.0, duhdx = 0.0, duhdy = 0.0;

      switch(_dim) {
        case 1:
          jac = dxdr[0];
          duhdx = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k);
          duhdx /= jac;
          break;
        case 2:
          jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];
          double drdx = dxds[1] / jac;
          double drdy = -dxds[0] / jac;
          double dsdx = -dxdr[1] / jac;
          double dsdy = dxdr[0] / jac;
          duhdx = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k) * drdx +
                  _intSpace->interpolateFieldAtQuadNode_sDerivative(_solution, k) * dsdx;
          duhdy = _intSpace->interpolateFieldAtQuadNode_rDerivative(_solution, k) * drdy +
                  _intSpace->interpolateFieldAtQuadNode_sDerivative(_solution, k) * dsdy;
          break;
      }

      norm[2] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - duhdx, 2);
      norm[3] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - duhdy, 2);
      norm[4] += J[nQuad * iElm + k] * w[k] *
                 (pow(dudxReconstruit - duhdx, 2) + pow(dudyReconstruit - duhdy, 2));
      if(solRefGrad) {
        std::vector<double> res(2, 0.);
        solRefGrad->eval(0, x, res);
        norm[5] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - res[0], 2);
        norm[6] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - res[1], 2);
        norm[7] += J[nQuad * iElm + k] * w[k] *
                   (pow(dudxReconstruit - res[0], 2) + pow(dudyReconstruit - res[1], 2));
        norm[18] += J[nQuad * iElm + k] * w[k] * pow(duhdx - res[0], 2);
        norm[19] += J[nQuad * iElm + k] * w[k] * pow(duhdy - res[1], 2);
        norm[20] += J[nQuad * iElm + k] * w[k] * (pow(duhdx - res[0], 2) + pow(duhdy - res[1], 2));
        fprintf(f, "%+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e\n", x[0], x[1],
                dudxReconstruit, res[0]);
      }
    }
  }
  norm[2] = sqrt(norm[2]);
  norm[3] = sqrt(norm[3]);
  norm[4] = sqrt(norm[4]);
  norm[5] = sqrt(norm[5]);
  norm[6] = sqrt(norm[6]);
  norm[7] = sqrt(norm[7]);
  norm[18] = sqrt(norm[18]);
  norm[19] = sqrt(norm[19]);
  norm[20] = sqrt(norm[20]);
  fclose(f);
}

void feRecovery::estimateHessError(std::vector<double> &norm, feVectorFunction *solRefHess)
{
  norm[8] = norm[9] = norm[10] = norm[11] = norm[12] = norm[13] = 0.0;

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimDerivation, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  std::vector<double> &solVec = _sol->getSolutionReference();

  // for(int iElm = 0; iElm < _nElm; ++iElm) {

  // std::set<int> &elemPatch = _patch->getPatch(4);

  FILE *f = fopen("d2Reconstruite.txt", "w");

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    // for(auto iElm : elemPatch) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }
    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // derivee reconstruite au point d'intégration : interpolation des valeurs issues des sommets
      double d2udx2Reconstruit = 0.0;
      double d2udxyReconstruit = 0.0;
      double d2udyxReconstruit = 0.0;
      double d2udy2Reconstruit = 0.0;
      // for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
      //   int v = _cnc->getNodeConnectivity(iElm, iNode);

      //   xLoc[0] = x[0] - _mesh->getVertex(v)->x();
      //   xLoc[1] = x[1] - _mesh->getVertex(v)->y();
      //   xLoc[2] = x[2] - _mesh->getVertex(v)->z();

      //   for(int i = 0; i < _dimDerivation; ++i) {
      //     monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
      //     d2udx2Reconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) *
      //     derivativeCoeff[v][0][i] * monomials[i]; d2udy2Reconstruit +=
      //     _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][1][i] * monomials[i];
      //   }
      // }

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // std::vector<double> &dux = derivativeCoeff[vNode][0];
        // std::vector<double> &duy = derivativeCoeff[vNode][1];
        // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
        // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
        // derivativeCoeff[vNode][0][i] * monomials[i]; d2udxyReconstruit +=
        // _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][i] * monomials[i];
        // d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
        // derivativeCoeff[vNode][2][i] * monomials[i]; d2udy2Reconstruit +=
        // _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][3][i] * monomials[i];
        // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i])
        // * pow(xLoc[1], _expY[i]);
        // }
        d2udx2Reconstruit +=
          _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][0];
        if(_dim == 2)
          d2udxyReconstruit +=
            _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][0];
        if(_dim == 2)
          d2udyxReconstruit +=
            _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][2][0];
        if(_dim == 2)
          d2udy2Reconstruit +=
            _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][3][0];
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
        // int v1, v2;
        // if(iEdge == _nEdgePerElm-1){
        //   v1 = _cnc->getNodeConnectivity(iElm, iEdge);
        //   v2 = _cnc->getNodeConnectivity(iElm, 0);
        // } else{
        //   v1 = _cnc->getNodeConnectivity(iElm, iEdge);
        //   v2 = _cnc->getNodeConnectivity(iElm, iEdge+1);
        // }
        // xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
        // xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
        // xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
        //   // Attention : iEdge+3 hardcoded for P2 interpolant
        //   // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
        //   _expX[i]) * pow(xLoc[1], _expY[i]);
        //   // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0],
        //   _expX[i]) * pow(xLoc[1], _expY[i]); monomials[i] = pow(xLoc[0], _expX[i]) *
        //   pow(xLoc[1], _expY[i]); d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3,
        //   k) * derivativeCoeffOnEdges[localEdge][0][0][i] * monomials[i]; d2udxyReconstruit +=
        //   _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        //   derivativeCoeffOnEdges[localEdge][0][1][i] * monomials[i]; d2udyxReconstruit +=
        //   _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        //   derivativeCoeffOnEdges[localEdge][0][2][i] * monomials[i]; d2udy2Reconstruit +=
        //   _intSpace->getFunctionAtQuadNode(iEdge+3, k) *
        //   derivativeCoeffOnEdges[localEdge][0][3][i] * monomials[i];
        // }
        if(_dim == 1)
          d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 2, k) *
                               derivativeCoeffOnEdges[localEdge][0][0][0];

        if(_dim == 2)
          d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                               derivativeCoeffOnEdges[localEdge][0][0][0];
        if(_dim == 2)
          d2udxyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                               derivativeCoeffOnEdges[localEdge][0][1][0];
        if(_dim == 2)
          d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                               derivativeCoeffOnEdges[localEdge][0][2][0];
        if(_dim == 2)
          d2udy2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) *
                               derivativeCoeffOnEdges[localEdge][0][3][0];
      }

      // std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      // std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      // _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      // _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      // double jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      // double drdx = dxds[1] / jac;
      // double drdy = -dxds[0] / jac;
      // double dsdx = -dxdr[1] / jac;
      // double dsdy = dxdr[0] / jac;

      // double duhdx = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdx +
      //                _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      // double duhdy = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdy +
      //                _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
      // norm[2] += J[nQuad * iElm + k] * w[k] * pow(d2udx2Reconstruit - duhdx, 2);
      // norm[3] += J[nQuad * iElm + k] * w[k] * pow(d2udy2Reconstruit - duhdy, 2);
      // norm[4] += J[nQuad * iElm + k] * w[k] * (pow(d2udx2Reconstruit - duhdx, 2) +
      // pow(d2udy2Reconstruit - duhdy, 2));

      if(solRefHess) {
        std::vector<double> res(4, 0.);
        solRefHess->eval(0, x, res);
        norm[8] += J[nQuad * iElm + k] * w[k] * pow(d2udx2Reconstruit - res[0], 2);
        norm[9] += J[nQuad * iElm + k] * w[k] * pow(d2udxyReconstruit - res[1], 2);
        norm[10] += J[nQuad * iElm + k] * w[k] *
                    (pow(d2udx2Reconstruit - res[0], 2) + pow(d2udxyReconstruit - res[1], 2));
        norm[11] += J[nQuad * iElm + k] * w[k] * pow(d2udyxReconstruit - res[2], 2);
        norm[12] += J[nQuad * iElm + k] * w[k] * pow(d2udy2Reconstruit - res[3], 2);
        norm[13] += J[nQuad * iElm + k] * w[k] *
                    (pow(d2udyxReconstruit - res[2], 2) + pow(d2udy2Reconstruit - res[3], 2));
        fprintf(f, "%+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e\n", x[0], x[1],
                d2udx2Reconstruit, res[0]);
      }
    }
  }
  norm[8] = sqrt(norm[8]);
  norm[9] = sqrt(norm[9]);
  norm[10] = sqrt(norm[10]);
  norm[11] = sqrt(norm[11]);
  norm[12] = sqrt(norm[12]);
  norm[13] = sqrt(norm[13]);
  fclose(f);
}

void feRecovery::estimated3Error(std::vector<double> &norm, feFunction *fund3udx)
{
  norm[21] = 0.0;

  std::vector<double> geoCoord(9, 0.);
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimDerivation, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  std::vector<double> &solVec = _sol->getSolutionReference();

  FILE *f = fopen("d3uReconstruite.txt", "w");

  // std::set<int> &elemPatch = _patch->getPatch(4);

  // for(int iElm = 0; iElm < _nElm; ++iElm) {
  for(int iElm = 0; iElm < _nElm - 3; ++iElm) {
    // for(auto v : _patch->getVertices()){
    // if(v < 8){
    // std::set<int> &elemPatch = _patch->getPatch(v);
    // norm[21] = 0.0;
    // for(auto iElm : elemPatch) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }

    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // derivee reconstruite au point d'intégration : interpolation des valeurs issues des sommets
      double d3udx3Reconstruit = 0.0;

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        d3udx3Reconstruit +=
          _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][0];
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        if(_dim == 1)
          d3udx3Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge + 2, k) *
                               derivativeCoeffOnEdges[localEdge][0][0][0];
      }

      if(fund3udx) {
        norm[21] += J[nQuad * iElm + k] * w[k] * pow(d3udx3Reconstruit - fund3udx->eval(0, x), 2);
      }

      fprintf(f, "%+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e\n", x[0], x[1],
              d3udx3Reconstruit, fund3udx->eval(0, x));
    }
  }

  int elmRef = _nElm - 3;

  // Last N elements : interpolate based on previous element
  for(int iElm = _nElm - 3; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm,
                                          _adr);
    for(size_t i = 0; i < _adr.size(); ++i) {
      _solution[i] = solVec[_adr[i]];
    }

    _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, geoCoord);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      double d3udx3Reconstruit = 0.0;

      int v = _cnc->getNodeConnectivity(elmRef, 1);

      xLoc[0] = x[0] - _mesh->getVertex(v)->x();
      xLoc[1] = x[1] - _mesh->getVertex(v)->y();
      xLoc[2] = x[2] - _mesh->getVertex(v)->z();

      for(int i = 0; i < _dimDerivation; ++i) {
        monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        d3udx3Reconstruit += derivativeCoeff[v][0][i] * monomials[i];
      }

      // // Vertices
      // for(int iVert = 0; iVert < _nVertPerElm; ++iVert){
      //   int vNode = _cnc->getNodeConnectivity(iElm, iVert);
      //   d3udx3Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) *
      //   derivativeCoeff[vNode][0][0];
      // }
      // // Edges
      // for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
      //   int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
      //   if(_dim == 1) d3udx3Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+2, k) *
      //   derivativeCoeffOnEdges[localEdge][0][0][0];
      // }

      if(fund3udx) {
        norm[21] += J[nQuad * iElm + k] * w[k] * pow(d3udx3Reconstruit - fund3udx->eval(0, x), 2);
      }

      fprintf(f, "%+-12.12e \t %+-12.12e \t %+-12.12e \t %+-12.12e\n", x[0], x[1],
              d3udx3Reconstruit, fund3udx->eval(0, x));
    }
  }

  norm[21] = sqrt(norm[21]);
  // }
  // }

  fclose(f);
}

static double NODAL_VALUES[6];
static double UVW[3];
static double SHAPE[6];
static double RESULT;

// Interpolates the recovered solution or its derivatives at point x.
double feRecovery::evalDerivative(int indexDerivative, double *x)
{
  // For P2 intspace only : 6 nodal values
  int elm = -1;
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    printf(
      "In feRecovery::evalDerivative : Warning - Point (%f, %f) was not found in the mesh.\n",
      x[0], x[1]);
    return 0.0;
  } else {
    NODAL_VALUES[0] = derivAtVertices[indexDerivative][_cnc->getNodeConnectivity(elm, 0)];
    NODAL_VALUES[1] = derivAtVertices[indexDerivative][_cnc->getNodeConnectivity(elm, 1)];
    NODAL_VALUES[2] = derivAtVertices[indexDerivative][_cnc->getNodeConnectivity(elm, 2)];
    NODAL_VALUES[3] = derivAtEdges[indexDerivative][fabs(_cnc->getEdgeConnectivity(elm, 0)) - 1];
    NODAL_VALUES[4] = derivAtEdges[indexDerivative][fabs(_cnc->getEdgeConnectivity(elm, 1)) - 1];
    NODAL_VALUES[5] = derivAtEdges[indexDerivative][fabs(_cnc->getEdgeConnectivity(elm, 2)) - 1];
    _intSpace->interpolateField(NODAL_VALUES, 6, UVW, SHAPE, RESULT);
    return RESULT;
  }

  // std::string meshName = "foo.msh";
  // bool curved = false;
  // bool verbose = false;
  // feMesh2DP1 *nmesh = new feMesh2DP1(meshName, curved, feMesh2DP1::mapType(), verbose);

  // // Sanity check :
  // FILE *f = fopen("test.pos","w");
  // fprintf(f, "View \"test\"{\n");

  // for(auto *t : nmesh->_elements){
  //   std::vector<double> val(3,0.0);
  //   std::vector<double> xpos(3,0.0);
  //   std::vector<double> ypos(3,0.0);

  //   for(int i = 0; i < 3; ++i){
  //     Vertex *v = t->getVertex(i);
  //     std::vector<double> x = {v->x(), v->y(), 0.0};
  //     std::vector<double> u(3,0.0);
  //     int elm;
  //     static_cast<feMesh2DP1*>(_mesh)->locateVertex(x,elm,u);

  //     std::vector<double> nodalValues(6,0.0);
  //     nodalValues[0] = derivAtVertices[0][_cnc->getNodeConnectivity(elm, 0)];
  //     nodalValues[1] = derivAtVertices[0][_cnc->getNodeConnectivity(elm, 1)];
  //     nodalValues[2] = derivAtVertices[0][_cnc->getNodeConnectivity(elm, 2)];
  //     nodalValues[3] = derivAtEdges[0][fabs(_cnc->getEdgeConnectivity(elm, 0))-1];
  //     nodalValues[4] = derivAtEdges[0][fabs(_cnc->getEdgeConnectivity(elm, 1))-1];
  //     nodalValues[5] = derivAtEdges[0][fabs(_cnc->getEdgeConnectivity(elm, 2))-1];

  //     val[i] = _intSpace->interpolateField(nodalValues,u.data());
  //     xpos[i] = v->x();
  //     ypos[i] = v->y();
  //   }

  //   fprintf(f,"ST(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%.16g, %.16g,
  //   %.16g};\n",
  //     xpos[0], ypos[0], 0., xpos[1], ypos[1], 0., xpos[2], ypos[2], 0.,
  //     val[0], val[1], val[2]);
  // }

  // fprintf(f, "};");
  // fclose(f);
}

// Dump all the recovered coefficients into a file to reload later
void feRecovery::writeRecovery(std::string fileName){
  FILE *f = fopen(fileName.c_str(), "w");
  if(f != nullptr) {

    fprintf(f, "%ld\n", derivAtVertices.size()); // Number of recovered derivatives at vertices

    for(size_t i = 0; i < derivAtVertices.size(); ++i){
      feInfo("Printing recovery %ld/%ld on %d DOFs", i+1, derivAtVertices.size(), derivAtVertices[i].size());
      fprintf(f, "%ld\n", derivAtVertices[i].size()); // Number of DOFs
      for(int j = 0; j < derivAtVertices[i].size(); ++j){
        fprintf(f, "%+-1.17e\n", derivAtVertices[i][j]);
      }
    }

    fprintf(f, "%ld\n", derivAtEdges.size()); // Number of recovered derivatives on the edges
    for(size_t i = 0; i < derivAtEdges.size(); ++i){
      feInfo("Printing recovery %ld/%ld on %ld DOFs", i+1, derivAtEdges.size(), derivAtEdges[i].size());
      fprintf(f, "%ld\n", derivAtEdges[i].size()); // Number of DOFs
      for(int j = 0; j < derivAtEdges[i].size(); ++j){
        fprintf(f, "%+-1.17e\n", derivAtEdges[i][j]);
      }
    }

    fclose(f);
  }
}
