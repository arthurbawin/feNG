#include "feRecovery.h"

#include <algorithm>

static bool firstRecovery = true;

// The names of the reconstructed fields
static std::map<std::pair<int, int>, std::string> suffix = {
  {{0, 0}, ""},     {{1, 0}, "dx"},   {{1, 1}, "dy"},   {{2, 0}, "dxx"},  {{2, 1}, "dxy"},
  {{2, 2}, "dyx"},  {{2, 3}, "dyy"},  {{3, 0}, "dxxx"}, {{3, 1}, "dxxy"}, {{3, 2}, "dxyx"},
  {{3, 3}, "dxyy"}, {{3, 4}, "dyxx"}, {{3, 5}, "dyxy"}, {{3, 6}, "dyyx"}, {{3, 7}, "dyyy"},
};

fePatch::fePatch(feCncGeo *cnc) {
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getNodeConnectivityCopy();
  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  _nVertices = _vertices.size();
  _nNodePerElm = cnc->getNbNodePerElem();
  std::vector<int> &connecNodes = cnc->getNodeConnectivityRef();

  int nElm = cnc->getNbElm();
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < _nNodePerElm; ++j) {
      vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i); // TODO : A verifier
    }
  }

  // Increase the size of the patch if there are too few elements
  for(auto &v : vertToElems){
    if(v.second.size() <= 2){
      std::vector<int> toAdd(0);
      for(auto e : v.second){
        for(int j = 0; j < _nNodePerElm; ++j) {
          for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]){
            toAdd.push_back(e2);
            // v.second.insert(e2);
          }
        }
      }
      for(auto e : toAdd)
        v.second.insert(e);
    }
  }
}

feRecovery::feRecovery(feMetaNumber *metaNumber, feSpace *space, feMesh *mesh, feSolution *sol,
                       std::vector<double> &norm, feFunction *solRef, std::string meshName, std::string metricMeshName,
                       feVectorFunction *solRefGrad)
  : _metaNumber(metaNumber), _mesh(mesh), _sol(sol), _intSpace(space) {
  _cnc = space->getCncGeo();
  _nElm = _cnc->getNbElm();
  _nNodePerElm = _cnc->getNbNodePerElem();
  _geoSpace = _cnc->getFeSpace();
  _degSol = space->getPolynomialDegree();
  // _degRec = _degSol + 1;
  _patch = new fePatch(_cnc);

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
  if(_dim == 2 && _cnc->getForme() == "TriP1") {
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

  allocateStructures();

  // To export all derivatives to a file, so that we don't need to store them all
  // std::string derivativesFileName = "derivatives.msh";
  // FILE *dFile = fopen(derivativesFileName.c_str(), "w");

  std::cout<<"Derivatives will be written in "<< metricMeshName <<std::endl;
  std::filebuf fbIn, fbOut;
  fbIn.open(meshName, std::ios::in);
  std::istream input(&fbIn);
  // if(firstRecovery){
    fbOut.open(metricMeshName, std::ios::out);
  // } else{
  //   // Append the derivatives of the next recoveries
  //   fbOut.open(metricMeshName, std::ios::app);
  // }

  std::ostream output(&fbOut);
  
  // if(firstRecovery){
    // firstRecovery = false;
    std::string buffer;
    // Copy .msh file
    while(getline(input, buffer)) {
      output << buffer << std::endl;
    }
    fbIn.close();
  // }

  for(int iDerivative = 0; iDerivative < _degSol + 1; ++iDerivative) {
    bool recoverDerivative = (iDerivative > 0);

    for(int i = 0; i < pow(_dim, iDerivative); ++i) {
      // Recovery of the solution if iDerivative = 0, of the derivatives if > 0
      solveLeastSquare(i, recoverDerivative);
      std::cout<<"Solving LEAST SQUARES"<<std::endl;
    }

    for(int i = 0; i < pow(_dim, iDerivative); ++i) {
      derivative(i, iDerivative, output);
      printf("Computed derivatives of solution of order %d\n", iDerivative);
    }

    // if(iDerivative == 0){
    //   estimateError(norm, solRef); // Compute L2 norm after the first reconstruction (the one for u)
    // } else if(iDerivative == 1){
    //   estimateH1Error(norm, solRefGrad);
    // }
    if(iDerivative == 0 && solRefGrad != nullptr){
      estimateH1Error(norm, solRefGrad);
    }
  }

  fbOut.close();

  getErrorPolynomials();

  freeStructures();
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

void feRecovery::freeStructures() {
  // Free structures
  free(valb);
  free(valA);
  free(indJ);
  free(indI);
  PetscErrorCode ierr;
  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&res);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&c);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&b);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatDestroy(&A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void feRecovery::solveLeastSquare(int indRecovery, bool recoverDerivative) {
  PetscErrorCode ierr;

  if(_geoSpace->getQuadratureWeights() != _intSpace->getQuadratureWeights()) {
    printf("Error : mismatch in number of quadrature points\n");
  }
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.0);
  std::vector<double> xLoc(3, 0.0);
  std::vector<double> monomials(_dimRecovery, 0.);

  // La solution interpolée aux points d'intégration
  double u, jac;

  int nQuad = _geoSpace->getNbQuadPoints();

  std::vector<int> &vertices = _patch->getVertices();
  for(auto v : vertices) {
    // Reset matrix and RHS
    ierr = MatZeroEntries(A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = VecZeroEntries(b);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // Get patch of elements
    std::set<int> &elemPatch = _patch->getPatch(v);

    // std::string name = "patchesVert" + std::to_string(v) + ".pos";

    // FILE *file = fopen(name.c_str(),"w");
    // std::string text = "View\"patchVertex" + std::to_string(v) + "\"{\n";
    // fprintf(file, text.c_str());

    for(auto elem : elemPatch) {
      _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()),
                                            elem);
      _intSpace->initializeSolution(_sol);
      geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), elem);

      // Sommets de l'element
      double r1[3] = {0., 0., 0.};
      double r2[3] = {1., 0., 0.};
      double r3[3] = {0., 1., 0.};
      std::vector<double> xv1(3, 0.0);
      std::vector<double> xv2(3, 0.0);
      std::vector<double> xv3(3, 0.0);
      _geoSpace->interpolateVectorField(geoCoord, r1, xv1);
      _geoSpace->interpolateVectorField(geoCoord, r2, xv2);
      _geoSpace->interpolateVectorField(geoCoord, r3, xv3);
       // fprintf(file,
       //    "ST(%g,%g,%g,%g,%g,%g,%g,%g,%g){%g,%g,%g};\n",
       //    xv1[0], xv1[1], 0.0, xv2[0], xv2[1], 0.0, xv3[0], xv3[1], 0.0, elem, elem, elem);

      // Loop over quad points and increment least square matrix
      for(int k = 0; k < nQuad; ++k) {
        jac = J[nQuad * elem + k];

        // TODO : normaliser ?
        _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
        xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        xLoc[2] = x[2] - _mesh->getVertex(v)->z();

        for(int i = 0; i < _dimRecovery; ++i)
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);

        for(int i = 0; i < _dimRecovery; ++i) {
          // valb[i] = jac * w[k] * u * monomials[i];
          for(int j = 0; j < _dimRecovery; ++j) {
            valA[_dimRecovery * i + j] = jac * w[k] * monomials[i] * monomials[j];
          }
        }
        ierr = MatSetValues(A, _dimRecovery, indI, _dimRecovery, indJ, valA, ADD_VALUES);

        if(recoverDerivative) {
          // The contributions to the derivative from the vertices must be evaluated and averaged at
          // quad nodes to avoid a trivial solution
          u = 0.0;
          for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
            int vNode = _cnc->getNodeConnectivity(elem, iNode);
            // Get the coefficients of the derivative (used only if recovering a derivative)
            std::vector<double> &du = derivativeCoeff[vNode][indRecovery];
            // printf("iNode = %d \t vertex = %d\n", iNode, vNode);
            xLoc[0] = x[0] - _mesh->getVertex(vNode)->x();
            xLoc[1] = x[1] - _mesh->getVertex(vNode)->y();
            xLoc[2] = x[2] - _mesh->getVertex(vNode)->z();
            // printf("iVert = %d - k = %2d - iNode = %d - xInt = %+-4.4f yInt = %+-4.4f xLoc =
            // %+-4.4f yLoc = %+-4.4f\n", v, k, iNode, x[0], x[1], xLoc[0], xLoc[1]);
            for(int i = 0; i < _dimDerivation; ++i) {
              u += _geoSpace->getFunctionAtQuadNode(iNode, k) * du[i] * pow(xLoc[0], _expX[i]) *
                   pow(xLoc[1], _expY[i]);
              // printf("Rec %d : %+-4.4f \t %+-4.4f \t %+-4.4f \t %+-4.4f\n", indRecovery, u,
              // _geoSpace->getFunctionAtQuadNode(iNode,k), du[i], pow(xLoc[0],
              // _expX[i])*pow(xLoc[1], _expY[i]));
            }
          }
        } else {
          // Simply interpolate the solution at quad nodes
          u = _intSpace->interpolateSolutionAtQuadNode(k);
        }

        for(int i = 0; i < _dimRecovery; ++i) valb[i] = jac * w[k] * u * monomials[i];

        ierr = VecSetValues(b, _dimRecovery, indI, valb, ADD_VALUES);
      }
    }

    // fprintf(file, "};\n");
    // fclose(file);

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);

    // Solve system
    PetscInt its;
    ierr = KSPSolve(ksp, b, c);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    KSPGetIterationNumber(ksp, &its);
    // PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
    VecSet(res, 0.0);
    MatMult(A, c, res);
    VecAXPY(res, -1.0, b);
    double normeAxmb = 0.0;
    ierr = VecNorm(res, NORM_MAX, &normeAxmb);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // std::cout<<"Norme du résidu matriciel Ax-b : "<<normeAxmb<<std::endl;

    // printf("Solution at vertex %f - %f - %f : \n", _mesh->getVertex(v)->x(),
    // _mesh->getVertex(v)->y(), _mesh->getVertex(v)->z()); VecView(c,PETSC_VIEWER_STDOUT_WORLD);

    // Copy coefficients in main structure
    PetscScalar *array;
    VecGetArray(c, &array);
    std::vector<double> coeffs((double *)array, (double *)array + _dimRecovery);
    VecRestoreArray(c, &array);
    recoveryCoeff[v][indRecovery] = coeffs;
  }
}

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

      // printf("Derivee à l'emplacement %d et %d = %+-4.4f - %+-4.4f\n", 2*indRecovery, 2*indRecovery+1, derivativeCoeff[v][2 * indRecovery][0], derivativeCoeff[v][2 * indRecovery+1][0]);
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
    int cnt = 1;
    for(auto v : vertices) {
      output << cnt++ << " " << recoveryCoeff[v][indRecovery][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }
  // Also write the derivatives to the common file
  for(int i = 0; i < _dim; ++i) {
    output << "$NodeData\n";
    std::string fieldName = "d" + std::to_string(iDerivative + 1) + _intSpace->getFieldID() +
                            suffix[{iDerivative + 1, _dim * indRecovery + i}];
    output << "1\n\"" << fieldName << "\"\n1\n3000\n3\n0\n1\n" << vertices.size() << "\n";
    int cnt = 1;
    for(auto v : vertices) {
      output << cnt++ << " " << derivativeCoeff[v][_dim * indRecovery + i][0] << std::endl;
    }
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

  for(int iElm = 0; iElm < _nElm; ++iElm) {
    _intSpace->initializeAddressingVector(_metaNumber->getNumbering(_intSpace->getFieldID()), iElm);
    _intSpace->initializeSolution(_sol);
    geoCoord = _mesh->getCoord(_intSpace->getCncGeoTag(), iElm);

    for(int k = 0; k < nQuad; ++k) {
      // Coordonnées des points d'intégration
      _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);

      // u reconstruit au point d'intégration : interpolation des valeurs issues des sommets
      double uReconstruit = 0.0;
      for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
        int v = _cnc->getNodeConnectivity(iElm, iNode);

        xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        xLoc[2] = x[2] - _mesh->getVertex(v)->z();

        for(int i = 0; i < _dimRecovery; ++i) {
          monomials[i] = pow(xLoc[0], _expXRecovery[i]) * pow(xLoc[1], _expYRecovery[i]);
          uReconstruit +=
            _geoSpace->getFunctionAtQuadNode(iNode, k) * recoveryCoeff[v][0][i] * monomials[i];
        }
      }

      norm[0] += J[nQuad * iElm + k] * w[k] *
                 pow(uReconstruit - _intSpace->interpolateSolutionAtQuadNode(k), 2);
      if(solRef) norm[1] += J[nQuad * iElm + k] * w[k] * pow(uReconstruit - solRef->eval(0, x), 2);
    }
  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
}

void feRecovery::estimateH1Error(std::vector<double> &norm, feVectorFunction *solRefGrad) {
  norm[0] = norm[1] = norm[2] = norm[3] = 0.0;

  std::vector<double> geoCoord;
  std::vector<double> x(3, 0.);
  std::vector<double> xLoc(3, 0.);
  std::vector<double> monomials(_dimRecovery, 0.);

  int nQuad = _geoSpace->getNbQuadPoints();
  std::vector<double> &w = _geoSpace->getQuadratureWeights();
  std::vector<double> &J = _cnc->getJacobians();

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
      for(int iNode = 0; iNode < _nNodePerElm; ++iNode) {
        int v = _cnc->getNodeConnectivity(iElm, iNode);

        xLoc[0] = x[0] - _mesh->getVertex(v)->x();
        xLoc[1] = x[1] - _mesh->getVertex(v)->y();
        xLoc[2] = x[2] - _mesh->getVertex(v)->z();

        for(int i = 0; i < _dimDerivation; ++i) {
          monomials[i] = pow(xLoc[0], _expX[i]) * pow(xLoc[1], _expY[i]);
          dudxReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][0][i] * monomials[i];
          dudyReconstruit += _geoSpace->getFunctionAtQuadNode(iNode, k) * derivativeCoeff[v][1][i] * monomials[i];
        }
      }

      norm[0] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - _intSpace->interpolateSolutionAtQuadNode_rDerivative(k), 2);
      norm[1] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - _intSpace->interpolateSolutionAtQuadNode_sDerivative(k), 2);
      if(solRefGrad){
        std::vector<double> res(2,0.);
        solRefGrad->eval(0,x,res);
        norm[2] += J[nQuad * iElm + k] * w[k] * pow(dudxReconstruit - res[0], 2);
        norm[3] += J[nQuad * iElm + k] * w[k] * pow(dudyReconstruit - res[1], 2);
      }
    }
  }
  norm[0] = sqrt(norm[0]);
  norm[1] = sqrt(norm[1]);
  norm[2] = sqrt(norm[2]);
  norm[3] = sqrt(norm[3]);
}