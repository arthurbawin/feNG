#include "feRecovery.h"

#include <algorithm>

#if defined(HAVE_EIGEN)
#include "../contrib/Eigen/Dense"
#endif

static inline double matNorm2(const std::vector<double> &v1, const std::vector<double> &v2, int n){
  double sqr = 0;
  for(int i = 0; i < n; i++) {
    sqr += (v1[i] - v2[i])*(v1[i] - v2[i]);
  }
  return sqrt(sqr);
}

// The names of the reconstructed fields
static std::map<std::pair<int, int>, std::string> suffix = {
  {{0, 0}, ""},     {{1, 0}, "dx"},   {{1, 1}, "dy"},   {{2, 0}, "dxx"},  {{2, 1}, "dxy"},
  {{2, 2}, "dyx"},  {{2, 3}, "dyy"},  {{3, 0}, "dxxx"}, {{3, 1}, "dxxy"}, {{3, 2}, "dxyx"},
  {{3, 3}, "dxyy"}, {{3, 4}, "dyxx"}, {{3, 5}, "dyxy"}, {{3, 6}, "dyyx"}, {{3, 7}, "dyyy"},
};

fePatch::fePatch(feCncGeo *cnc, feMesh *mesh) {
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getNodeConnectivityCopy();
  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  _nVertices = _vertices.size();
  _nNodePerElm = cnc->getNbNodePerElem();
  _nEdgePerElm = cnc->getNbEdgePerElem();
  std::vector<int> &connecNodes = cnc->getNodeConnectivityRef();
  std::vector<int> &connecEdges = cnc->getEdgeConnectivityRef();

  int nElm = cnc->getNbElm();
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < _nNodePerElm; ++j) {
      vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i); // TODO : A verifier
    }
  }

  for(auto e : mesh->_edges){
    // Insert the patches of both edge vertices to the edge's patch
    int v0 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(0));
    int v1 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(1));
    edgeToElems[e.getTag()].insert(vertToElems[v0].begin(), vertToElems[v0].end());
    edgeToElems[e.getTag()].insert(vertToElems[v1].begin(), vertToElems[v1].end());
  }

  // Increase the size of the patch if there are too few elements
  for(auto &v : vertToElems) {
    if(v.second.size() <= 2) {
      std::vector<int> toAdd(0);
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
}

feRecovery::feRecovery(feMetaNumber *metaNumber, feSpace *space, feMesh *mesh, feSolution *sol,
                       std::vector<double> &norm, feFunction *solRef, std::string meshName,
                       std::string metricMeshName, feVectorFunction *solRefGrad, feVectorFunction *solRefHess, bool append)
  : _metaNumber(metaNumber), _mesh(mesh), _sol(sol), _intSpace(space) {
  _cnc = space->getCncGeo();
  _nElm = _cnc->getNbElm();
  _nNodePerElm = _cnc->getNbNodePerElem();
  _geoSpace = _cnc->getFeSpace();
  _degSol = space->getPolynomialDegree();
  // _degRec = _degSol + 1;
  _patch = new fePatch(_cnc, _mesh);

  _dim = mesh->getDim();

  // Set total number of recoveries/derivations
  _nTotalRecoveries = pow(_dim, _degSol); // TODO : A verifier
  _nTotalDerivations = pow(_dim, _degSol + 1);

  // The dimension of the polynomial bases for recoveries is the one of degree k+1 :
  if(_dim == 1) {
    printf("Error : No recovery method available for one-dimensional mesh.\n");
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
      // } else if(space->getCncGeo()->getForme() == "QuadP1" || space->getCncGeo()->getForme() ==
      // "QuadP2"){
      //   _dimRecovery   = (_degSol+2)*(_degSol+2)*(_degSol+2);
      //   _dimDerivation = (_degSol+1)*(_degSol+1)*(_degSol+1);
    } else {
      printf("Error : Mesh connectivity \"%s\" not available for polynomial recovery.\n",
             _cnc->getForme().c_str());
    }
  }

  // The exponents of the monomials :
  _expXRecovery.resize(_dimRecovery, 0);
  _expYRecovery.resize(_dimRecovery, 0);
  _expZRecovery.resize(_dimRecovery, 0);
  _expX.resize(_dimDerivation, 0);
  _expY.resize(_dimDerivation, 0);
  _expZ.resize(_dimDerivation, 0);
  int ind = 0, n = _degSol + 1;
  if(_dim == 2 && (_cnc->getForme() == "TriP1" || _cnc->getForme() == "TriP2")) {
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
  } else {
    printf(
      "TODO : Error : Recovery coefficients not implemented for chosen dimension and/or element\n");
  }

  // allocateStructures();
  printf("Computing inverses with Eigen...\n");
  matrixInverseEigen();
  printf("Done\n");

  // To export all derivatives to a file, so that we don't need to store them all
  // std::string derivativesFileName = "derivatives.msh";
  // FILE *dFile = fopen(derivativesFileName.c_str(), "w");

  printf("Info in feRecovery : Derivatives will be written to file \"%s\"\n", metricMeshName.c_str());
  std::filebuf fbIn, fbOut;
  fbIn.open(meshName, std::ios::in);
  std::istream input(&fbIn);
  fbOut.open(metricMeshName, std::ios::out);
  std::ostream output(&fbOut);

  std::string buffer;
  if(append){
    // Copy .msh file
    while(getline(input, buffer)) { 
      output << buffer << std::endl;
    }
    fbIn.close();
  } else{
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

  for(int iDerivative = 0; iDerivative < _degSol + 1; ++iDerivative) {
    bool recoverDerivative = (iDerivative > 0);

    // for(int i = 0; i < pow(_dim, iDerivative); ++i) {
      // Recovery of the solution if iDerivative = 0, of the derivatives if > 0
      solveLeastSquareEigen(pow(_dim, iDerivative), recoverDerivative);
      // solveLeastSquare(i, recoverDerivative);
    // }

    for(int i = 0; i < pow(_dim, iDerivative); ++i) {
      derivative(i, iDerivative, output);
      printf("Computed derivatives of solution of order %d\n", iDerivative);
    }

    // if(iDerivative == 0){
    //   estimateError(norm, solRef); // Compute L2 norm after the first reconstruction (the one for
    //   u)
    // } else if(iDerivative == 1){
    //   estimateH1Error(norm, solRefGrad);
    // }
    if(iDerivative == 0 && solRef != nullptr) { estimateError(norm, solRef); }
    if(iDerivative == 0 && solRefGrad != nullptr) { estimateH1Error(norm, solRefGrad); }
    if(iDerivative == 1 && solRefGrad != nullptr) { estimateDudxError(norm, solRefGrad); }
    if(iDerivative == 1 && solRefHess != nullptr) { estimateHessError(norm, solRefHess); }
  }

  fbOut.close();

  getErrorPolynomials();

  // freeStructures();

  // std::string gmshStr = "gmsh " + metricMeshName + " &";
  // system(gmshStr.c_str());
}

void feRecovery::allocateStructures() {
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
  for(int i = 0; i < _dimRecovery; ++i) { indI[i] = indJ[i] = i; }
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

void feRecovery::matrixInverseEigen(){

  // if(_degSol != 2){
  //   printf("Modifier les matrices eigen pour deg != 2\n");
  //   exit(-1);
  // }

  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord, x(3, 0.0), xLoc(3, 0.0), monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  // Matrices defined on the vertices
  for(auto v : _patch->getVertices()) {

    Eigen::Matrix<double, 6, 6> myMat6 = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 10, 10> myMat10 = Eigen::Matrix<double, 10, 10>::Zero();
    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      // std::cout<<elem<<std::endl;
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem);
      _intSpace->initializeSolution(_sol);
      geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

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
        xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        xLoc[2] = x[2] - _mesh->getVertex(v)->z();

        for(int i = 0; i < _dimRecovery; ++i){
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
        }

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol){
              case 1: myMat6(i,j) += jac * w[k] * monomials[i] * monomials[j]; break;
              case 2: myMat10(i,j) += jac * w[k] * monomials[i] * monomials[j]; break;
            }
          }
        }
      }
    }
    switch(_degSol){
      case 1: lsInvAtVertices6[v] = myMat6.inverse(); break;
      case 2: lsInvAtVertices10[v] = myMat10.inverse(); break;
    }
  }

  // Matrices defined on the edges
  for(auto e : _mesh->_edges) {

    // TODO : boucler sur le nombre de DOFS par edge, ici on suppose juste un P2 avec 1 dof

    Eigen::Matrix<double, 6, 6> myMat6 = Eigen::Matrix<double, 6, 6>::Zero();
    Eigen::Matrix<double, 10, 10> myMat10 = Eigen::Matrix<double, 10, 10>::Zero();

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      // std::cout<<elem<<std::endl;
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem);
      _intSpace->initializeSolution(_sol);
      geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        double jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        // printf("%d : (%f,%f) et %d : (%f,%f)\n", 
        //   e.getTag(0), _mesh->getVertex(e.getTag(0))->x(), _mesh->getVertex(e.getTag(0))->y(),
        //   e.getTag(1), _mesh->getVertex(e.getTag(1))->x(), _mesh->getVertex(e.getTag(1))->y());
        xLoc[0] = x[0] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->x() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->x())/2.0;
        xLoc[1] = x[1] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->y() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->y())/2.0;
        xLoc[2] = x[2] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->z() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->z())/2.0;

        for(int i = 0; i < _dimRecovery; ++i)
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);

        for(int i = 0; i < _dimRecovery; ++i) {
          for(int j = 0; j < _dimRecovery; ++j) {
            switch(_degSol){
              case 1: myMat6(i,j) += jac * w[k] * monomials[i] * monomials[j]; break;
              case 2: myMat10(i,j) += jac * w[k] * monomials[i] * monomials[j]; break;
            }
          }
        }
      }
    }
    switch(_degSol){
      case 1: lsInvAtVertices6OnEdges[e.getTag()] = myMat6.inverse(); break;
      case 2: lsInvAtVertices10OnEdges[e.getTag()] = myMat10.inverse(); break;
    }
  }
}

void feRecovery::solveLeastSquareEigen(int indRecovery, bool recoverDerivative) {

  printf("Calling solveLS with %d\n", indRecovery);

  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.0);
  std::vector<double> xLoc(3, 0.0);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();

  std::vector<int> &vertices = _patch->getVertices();

  std::vector<double> u(indRecovery, 0.);

  int nDOFPerElem = _intSpace->getNbFunctions();
  int _nEdgePerElm = _cnc->getNbEdgePerElem();
  int _nVertPerElm = _nNodePerElm;
  
  for(auto v : vertices) {

    std::vector<Eigen::VectorXd> RHS6(indRecovery, Eigen::MatrixXd::Zero(6,1));
    std::vector<Eigen::VectorXd> RHS10(indRecovery, Eigen::MatrixXd::Zero(10,1));

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem);
      _intSpace->initializeSolution(_sol);
      geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

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
        xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        xLoc[2] = x[2] - _mesh->getVertex(v)->z();

        for(int i = 0; i < _dimRecovery; ++i){
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
        }

        if(recoverDerivative) {
          // The contributions to the derivative from the vertices must be evaluated and averaged at
          // quad nodes to avoid a trivial solution
          // u = 0.0;
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
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
          //       u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0], _expX[i]) *
          //            pow(xLoc[1], _expY[i]);
          //     }
          //   }
          // }

          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
            // Vertices
            // printf("elem with vertices\n");
            for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
              int vNode = _cnc->getNodeConnectivity(elem, iVert);
              std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
              // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
              // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
              // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
              // printf("%d en %f - %f\n", vNode, _mesh->getVertex(vNode)->x(), _mesh->getVertex(vNode)->y());
              // for(int i = 0; i < _dimDerivation; ++i) {
                // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              // }
            }
            // Edges
            for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
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
              //   v1, v2, _mesh->getVertex(v1)->x(), _mesh->getVertex(v1)->y(), _mesh->getVertex(v2)->x(), _mesh->getVertex(v2)->y());

              // Edge e(_mesh->getVertex(v1),_mesh->getVertex(v2));
              // std::set<Edge, EdgeLessThan>::iterator ret;
              // ret = _mesh->_edges.find(e);
              // if(ret != _mesh->_edges.end()){
              //   printf("sanityCheck : edge %d was found at %d : relie les sommets %d et %d\n", localEdge, ret->getTag(), ret->getTag(0), ret->getTag(1));
              // } else{
              //   printf("sanityCheck : edge %d was NOT found\n", localEdge, _mesh->_edges);
              // }
              // for(int i = 0; i < _dimDerivation; ++i) {
                // Attention : iEdge+3 hardcoded for P2 interpolant
                // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[0];
              // }
            }
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
            u[iDeriv] = _intSpace->interpolateSolutionAtQuadNode(k);
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
          for(int i = 0; i < _dimRecovery; ++i){
            switch(_degSol){
              case 1: RHS6[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i]; break; 
              case 2: RHS10[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i]; break; 
            }
          }
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
      switch(_degSol){
        case 1: 
        {
          Eigen::VectorXd sol = lsInvAtVertices6[v]*RHS6[iDeriv];
          recoveryCoeff[v][iDeriv] = std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2:
        {
          Eigen::VectorXd sol = lsInvAtVertices10[v]*RHS10[iDeriv];
          recoveryCoeff[v][iDeriv] = std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
      }
      
    }
  }

  for(auto e : _mesh->_edges) {

    std::vector<Eigen::VectorXd> RHS6(indRecovery, Eigen::MatrixXd::Zero(6,1));
    std::vector<Eigen::VectorXd> RHS10(indRecovery, Eigen::MatrixXd::Zero(10,1));

    // Get patch of elements
    std::set<int> &elemPatch = _patch->getEdgePatch(e.getTag());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem);
      _intSpace->initializeSolution(_sol);
      geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

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
        xLoc[0] = x[0] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->x() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->x())/2.0;
        xLoc[1] = x[1] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->y() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->y())/2.0;
        xLoc[2] = x[2] - (_mesh->getVertexFromGmshNodeTag(e.getTag(0))->z() + _mesh->getVertexFromGmshNodeTag(e.getTag(1))->z())/2.0;

        for(int i = 0; i < _dimRecovery; ++i){
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
        }

        if(recoverDerivative) {
          // The contributions to the derivative from the vertices must be evaluated and averaged at
          // quad nodes to avoid a trivial solution
          // u = 0.0;
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
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
          //       u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          //     }
          //   }
          // }
          
          // Get the coefficients of the derivative (used only if recovering a derivative)
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
            // Vertices
            // printf("elem with vertices\n");
            for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
              int vNode = _cnc->getNodeConnectivity(elem, iVert);
              std::vector<double> &du = derivativeCoeff[vNode][iDeriv];
              // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
              // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
              // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
              // printf("%d en %f - %f\n", vNode, _mesh->getVertex(vNode)->x(), _mesh->getVertex(vNode)->y());
              // for(int i = 0; i < _dimDerivation; ++i) {
                // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              // }
            }
            // Edges
            for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
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
              //   v1, v2, _mesh->getVertex(v1)->x(), _mesh->getVertex(v1)->y(), _mesh->getVertex(v2)->x(), _mesh->getVertex(v2)->y());

              // Edge e(_mesh->getVertex(v1),_mesh->getVertex(v2));
              // std::set<Edge, EdgeLessThan>::iterator ret;
              // ret = _mesh->_edges.find(e);
              // if(ret != _mesh->_edges.end()){
              //   printf("sanityCheck : edge %d was found at %d : relie les sommets %d et %d\n", localEdge, ret->getTag(), ret->getTag(0), ret->getTag(1));
              // } else{
              //   printf("sanityCheck : edge %d was NOT found\n", localEdge, _mesh->_edges);
              // }
              // for(int i = 0; i < _dimDerivation; ++i) {
                // Attention : iEdge+3 hardcoded for P2 interpolant
                // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
                u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[0];
              // }
            }
          }

        } else {
          // Simply interpolate the solution at quad nodes
          for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
            u[iDeriv] = _intSpace->interpolateSolutionAtQuadNode(k);
          }
        }

        for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
          for(int i = 0; i < _dimRecovery; ++i){
            switch(_degSol){
              case 1: RHS6[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i]; break; 
              case 2: RHS10[iDeriv](i) += jac * w[k] * u[iDeriv] * monomials[i]; break; 
            }
          }
        }
      }
    }
    for(int iDeriv = 0; iDeriv < indRecovery; ++iDeriv){
      switch(_degSol){
        case 1: 
        {
          Eigen::VectorXd sol = lsInvAtVertices6OnEdges[e.getTag()]*RHS6[iDeriv];
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] = std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
          break;
        }
        case 2:
        {
          Eigen::VectorXd sol = lsInvAtVertices10OnEdges[e.getTag()]*RHS10[iDeriv];
          recoveryCoeffOnEdges[e.getTag()][0][iDeriv] = std::vector<double>(sol.data(), sol.data() + sol.rows() * sol.cols());
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
//                                             elem);
//       _intSpace->initializeSolution(_sol);
//       geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

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
//           // The contributions to the derivative from the vertices must be evaluated and averaged at
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

void feRecovery::derivative(int indRecovery, int iDerivative, std::ostream &output) {
  std::vector<int> &vertices = _patch->getVertices();
  if(_dim == 2) {
    for(auto v : vertices) {
      int indX = 0, indY = 0;
      std::vector<double> &u = recoveryCoeff[v][indRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      std::vector<double> dudy(_dimDerivation, 0.);
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) { dudx[indX++] = _expXRecovery[i] * u[i]; }
        if(_expYRecovery[i] != 0) { dudy[indY++] = _expYRecovery[i] * u[i]; }
      }
      derivativeCoeff[v][2 * indRecovery] = dudx;
      derivativeCoeff[v][2 * indRecovery + 1] = dudy;

      // printf("Derivee à l'emplacement %d et %d = %+-4.4f - %+-4.4f\n", 2*indRecovery,
      // 2*indRecovery+1, derivativeCoeff[v][2 * indRecovery][0], derivativeCoeff[v][2 *
      // indRecovery+1][0]);
    }

    for(auto e : _mesh->_edges) {
      int indX = 0, indY = 0;
      std::vector<double> &u = recoveryCoeffOnEdges[e.getTag()][0][indRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      std::vector<double> dudy(_dimDerivation, 0.);
      for(int i = 0; i < _dimRecovery; ++i) {
        if(_expXRecovery[i] != 0) { dudx[indX++] = _expXRecovery[i] * u[i]; }
        if(_expYRecovery[i] != 0) { dudy[indY++] = _expYRecovery[i] * u[i]; }
      }
      derivativeCoeffOnEdges[e.getTag()][0][2 * indRecovery] = dudx;
      derivativeCoeffOnEdges[e.getTag()][0][2 * indRecovery + 1] = dudy;
    }

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
    printf("TODO : Derivatives for recoveries in 0D, 1D or 3D\n");
  }

  // // Print the solution at first pass
  // if(iDerivative == 0){
  //   fprintf(dFile, "$NodeData\n");
  //   fprintf(dFile,"1\n\"%s\"\n1\n3000\n3\n0\n1\n%lu\n", _intSpace->getFieldID().c_str(),
  //   vertices.size()); int cnt = 1; for(auto v : vertices){
  //     fprintf(dFile, "%6u %12.5E\n", cnt++, recoveryCoeff[v][indRecovery][0]);
  //   }
  //   fprintf(dFile, "$EndNodeData\n");
  // }
  // // Also write the derivatives to the common file
  // for(int i = 0; i < _dim; ++i){
  //   fprintf(dFile, "$NodeData\n");
  //   std::string fieldName = "d" + std::to_string(iDerivative+1) + _intSpace->getFieldID() +
  //   suffix[{iDerivative+1, _dim*indRecovery+i}];
  //   fprintf(dFile,"1\n\"%s\"\n1\n3000\n3\n0\n1\n%lu\n", fieldName.c_str(), vertices.size());
  //   int cnt = 1;
  //   for(auto v : vertices){
  //     fprintf(dFile, "%6u %12.5E\n", cnt++, derivativeCoeff[v][_dim*indRecovery+i][0]);
  //   }
  //   fprintf(dFile, "$EndNodeData\n");
  // }

  // Print the solution at first pass
  if(iDerivative == 0) {
    output << "$NodeData\n";
    output << "1\n\"" << _intSpace->getFieldID() << "\"\n1\n3000\n3\n0\n1\n"
           << vertices.size() << "\n";
    int cnt = 0;
    for(auto v : vertices) {
      // printf("Writing %d at vertex %d\n", v, _mesh->getVertex(cnt)->getTag());
      // output << _mesh->getVertex(cnt++)->getTag() << " " << recoveryCoeff[v][indRecovery][0] << std::endl;
      output << _mesh->getVertex(v)->getTag() << " " << recoveryCoeff[v][indRecovery][0] << std::endl;
      // output << v << " " << recoveryCoeff[v][indRecovery][0] << std::endl;
    }
    // output << 67 << " " << 0.0 << std::endl;
    // output << 7624 << " " << 0.0 << std::endl;
    // output << 11500 << " " << 0.0 << std::endl;
    // output << 20240 << " " << 0.0 << std::endl;
    // output << 4536 << " " << 0.0 << std::endl;
    // output << 8419 << " " << 0.0 << std::endl;
    // output << 18654 << " " << 0.0 << std::endl;
    output << "$EndNodeData\n";
  }
  // Also write the derivatives to the common file
  for(int i = 0; i < _dim; ++i) {
    output << "$NodeData\n";
    std::string fieldName = "d" + std::to_string(iDerivative + 1) + _intSpace->getFieldID() +
                            suffix[{iDerivative + 1, _dim * indRecovery + i}];
    output << "1\n\"" << fieldName << "\"\n1\n3000\n3\n0\n1\n" << vertices.size() << "\n";
    int cnt = 0;
    for(auto v : vertices) {
      // output << _mesh->getVertex(cnt++)->getTag() << " " << derivativeCoeff[v][_dim * indRecovery + i][0] << std::endl;
      output << _mesh->getVertex(v)->getTag() << " " << derivativeCoeff[v][_dim * indRecovery + i][0] << std::endl;
    }
    // output << 67 << " " << 0.0 << std::endl;
    // output << 7624 << " " << 0.0 << std::endl;
    // output << 11500 << " " << 0.0 << std::endl;
    // output << 20240 << " " << 0.0 << std::endl;
    // output << 4536 << " " << 0.0 << std::endl;
    // output << 8419 << " " << 0.0 << std::endl;
    // output << 18654 << " " << 0.0 << std::endl;
    output << "$EndNodeData\n";
  }
}

void feRecovery::getErrorPolynomials() {
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
        error[1] = derivativeCoeff[v][1][0] + derivativeCoeff[v][2][0] + derivativeCoeff[v][4][0] +
                   derivativeCoeff[v][8][0];
        error[2] = derivativeCoeff[v][3][0] + derivativeCoeff[v][5][0] + derivativeCoeff[v][6][0] +
                   derivativeCoeff[v][9][0] + derivativeCoeff[v][10][0] + derivativeCoeff[v][12][0];
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

void feRecovery::estimateError(std::vector<double> &norm, feFunction *solRef) {
  norm[0] = norm[1] = 0.0;

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(_sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);

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
      //     uReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * recoveryCoeff[v][0][i] * monomials[i];
      //   }
      // }

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][i] * monomials[i];
          uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][0];
        // }
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        int v1, v2;
        if(iEdge == _nEdgePerElm-1){
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, 0);
        } else{
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, iEdge+1);
        }
        xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
        xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
        xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // uReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * recoveryCoeffOnEdges[localEdge][0][0][i] * monomials[i];
          uReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * recoveryCoeffOnEdges[localEdge][0][0][0];
        // }
      }

      norm[0] += J[nQuad * iElm + k] * w[k] * pow(uReconstruit - _intSpace->interpolateSolutionAtQuadNode(k), 2);
      if(solRef) norm[1] += J[nQuad * iElm + k] * w[k] * pow(uReconstruit - solRef->eval(0, x), 2);
    }
  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
}

void feRecovery::estimateDudxError(std::vector<double> &norm, feVectorFunction *solRefGrad) {
  norm[14] = norm[15] = norm[16] = norm[17] = 0.0;

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(_sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // u reconstruit au point d'intégration : interpolation des valeurs issues des sommets
      double dudxReconstruit = 0.0;
      double dudyReconstruit = 0.0;

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        // xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        // xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // uReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][i] * monomials[i];
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][0][0];
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * recoveryCoeff[vNode][1][0];
        // }
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
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
          // uReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * recoveryCoeffOnEdges[localEdge][0][0][i] * monomials[i];
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * recoveryCoeffOnEdges[localEdge][0][0][0];
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * recoveryCoeffOnEdges[localEdge][0][1][0];
        // }
      }

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      double jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      double drdx = dxds[1] / jac;
      double drdy = -dxds[0] / jac;
      double dsdx = -dxdr[1] / jac;
      double dsdy = dxdr[0] / jac;

      double duhdx = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdx + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      double duhdy = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdy + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;

      norm[14] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - duhdx, 2);
      norm[15] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - duhdy, 2);
      if(solRefGrad){
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

void feRecovery::estimateH1Error(std::vector<double> &norm, feVectorFunction *solRefGrad) {
  norm[2] = norm[3] = norm[4] = norm[5] = norm[6] = norm[7] = 0.0;
  norm[18] = norm[19] = norm[20] = 0.0;

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimDerivation, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(_sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);

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
      //     dudxReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][0][i] * monomials[i];
      //     dudyReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][1][i] * monomials[i];
      //   }
      // }

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // std::vector<double> &dux = derivativeCoeff[vNode][0];
        // std::vector<double> &duy = derivativeCoeff[vNode][1];
        xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // dudxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][i] * monomials[i];
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][0];
          // dudyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][i] * monomials[i];
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][0];
          // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // }
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
        int v1, v2;
        if(iEdge == _nEdgePerElm-1){
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, 0);
        } else{
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, iEdge+1);
        }
        xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
        xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
        xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
          // Attention : iEdge+3 hardcoded for P2 interpolant
          // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][0][i] * monomials[i];
          dudxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][0][0];
          // dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][1][i] * monomials[i];
          dudyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][1][0];
        // }
      }

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      double jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      double drdx = dxds[1] / jac;
      double drdy = -dxds[0] / jac;
      double dsdx = -dxdr[1] / jac;
      double dsdy = dxdr[0] / jac;

      double duhdx = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdx + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      double duhdy = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdy + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
      norm[2] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - duhdx, 2);
      norm[3] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - duhdy, 2);
      norm[4] += J[nQuad * iElm + k] * w[k] * (pow(dudxReconstruit - duhdx, 2) + pow(dudyReconstruit - duhdy, 2));
      if(solRefGrad) {
        std::vector<double> res(2, 0.);
        solRefGrad->eval(0, x, res);
        norm[5] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - res[0], 2);
        norm[6] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - res[1], 2);
        norm[7] += J[nQuad * iElm + k] * w[k] * (pow(dudxReconstruit - res[0], 2) + pow(dudyReconstruit - res[1], 2));
        norm[18] += J[nQuad * iElm + k] * w[k] * pow(duhdx - res[0], 2);
        norm[19] += J[nQuad * iElm + k] * w[k] * pow(duhdy - res[1], 2);
        norm[20] += J[nQuad * iElm + k] * w[k] * (pow(duhdx - res[0], 2) + pow(duhdy - res[1], 2));
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
}

void feRecovery::estimateHessError(std::vector<double> &norm, feVectorFunction *solRefHess) {
  norm[8] = norm[9] = norm[10] = norm[11] = norm[12] = norm[13] = 0.0;

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimDerivation, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  int _nVertPerElm = _nNodePerElm;
  int _nEdgePerElm = _cnc->getNbEdgePerElem();

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(_sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);

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
      //     d2udx2Reconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][0][i] * monomials[i];
      //     d2udy2Reconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][1][i] * monomials[i];
      //   }
      // }

      // Vertices
      for(int iVert = 0; iVert < _nVertPerElm; ++iVert){ 
        int vNode = _cnc->getNodeConnectivity(iElm, iVert);
        // std::vector<double> &dux = derivativeCoeff[vNode][0];
        // std::vector<double> &duy = derivativeCoeff[vNode][1];
        xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
        xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
        xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
        // for(int i = 0; i < _dimDerivation; ++i) {
          // monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          // d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][i] * monomials[i];
          // d2udxyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][i] * monomials[i];
          // d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][2][i] * monomials[i];
          // d2udy2Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][3][i] * monomials[i];
          // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        // }
        d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][0][0];
        d2udxyReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][1][0];
        d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][2][0];
        d2udy2Reconstruit += _intSpace->getFunctionAtQuadNode(iVert, k) * derivativeCoeff[vNode][3][0];
      }
      // Edges
      for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge){
        int localEdge = fabs(_cnc->getEdgeConnectivity(iElm, iEdge));
        // std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iDeriv];
        int v1, v2;
        if(iEdge == _nEdgePerElm-1){
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, 0);
        } else{
          v1 = _cnc->getNodeConnectivity(iElm, iEdge);
          v2 = _cnc->getNodeConnectivity(iElm, iEdge+1);
        }
        xLoc[0] = x[0] - (_mesh->getVertex(v1)->x() + _mesh->getVertex(v2)->x())/2.0;
        xLoc[1] = x[1] - (_mesh->getVertex(v1)->y() + _mesh->getVertex(v2)->y())/2.0;
        xLoc[2] = x[2] - (_mesh->getVertex(v1)->z() + _mesh->getVertex(v2)->z())/2.0;
        // for(int i = 0; i < _dimDerivation; ++i) {
        //   // Attention : iEdge+3 hardcoded for P2 interpolant
        //   // u[iDeriv] += _geoSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        //   // u[iDeriv] += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * du[i] * pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        //   monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
        //   d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][0][i] * monomials[i];
        //   d2udxyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][1][i] * monomials[i];
        //   d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][2][i] * monomials[i];
        //   d2udy2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][3][i] * monomials[i];
        // }
        d2udx2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][0][0];
        d2udxyReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][1][0];
        d2udyxReconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][2][0];
        d2udy2Reconstruit += _intSpace->getFunctionAtQuadNode(iEdge+3, k) * derivativeCoeffOnEdges[localEdge][0][3][0];
      }

      std::vector<double> dxdr(3, 0.0); // [dx/dr, dy/dr, dz/dr]
      std::vector<double> dxds(3, 0.0); // [dx/ds, dy/ds, dz/ds]
      _geoSpace->interpolateVectorFieldAtQuadNode_rDerivative(geoCoord, k, dxdr);
      _geoSpace->interpolateVectorFieldAtQuadNode_sDerivative(geoCoord, k, dxds);
      double jac = dxdr[0] * dxds[1] - dxdr[1] * dxds[0];

      double drdx = dxds[1] / jac;
      double drdy = -dxds[0] / jac;
      double dsdx = -dxdr[1] / jac;
      double dsdy = dxdr[0] / jac;

      double duhdx = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdx + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdx;
      double duhdy = _intSpace->interpolateSolutionAtQuadNode_rDerivative(k) * drdy + 
                     _intSpace->interpolateSolutionAtQuadNode_sDerivative(k) * dsdy;
      // norm[2] += J[nQuad * iElm + k] * w[k] * pow(d2udx2Reconstruit - duhdx, 2);
      // norm[3] += J[nQuad * iElm + k] * w[k] * pow(d2udy2Reconstruit - duhdy, 2);
      // norm[4] += J[nQuad * iElm + k] * w[k] * (pow(d2udx2Reconstruit - duhdx, 2) + pow(d2udy2Reconstruit - duhdy, 2));

      if(solRefHess) {
        std::vector<double> res(4, 0.);
        solRefHess->eval(0, x, res);
        norm[8]  += J[nQuad * iElm + k] * w[k] *  pow(d2udx2Reconstruit - res[0], 2);
        norm[9]  += J[nQuad * iElm + k] * w[k] *  pow(d2udxyReconstruit - res[1], 2);
        norm[10] += J[nQuad * iElm + k] * w[k] * (pow(d2udx2Reconstruit - res[0], 2) + pow(d2udxyReconstruit - res[1], 2));
        norm[11] += J[nQuad * iElm + k] * w[k] *  pow(d2udyxReconstruit - res[2], 2);
        norm[12] += J[nQuad * iElm + k] * w[k] *  pow(d2udy2Reconstruit - res[3], 2);
        norm[13] += J[nQuad * iElm + k] * w[k] * (pow(d2udyxReconstruit - res[2], 2) + pow(d2udy2Reconstruit - res[3], 2));
      }
    }
  }
  norm[8] = sqrt(norm[8]);
  norm[9] = sqrt(norm[9]);
  norm[10] = sqrt(norm[10]);
  norm[11] = sqrt(norm[11]);
  norm[12] = sqrt(norm[12]);
  norm[13] = sqrt(norm[13]);
}