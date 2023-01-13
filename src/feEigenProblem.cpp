#include "feEigenProblem.h"
#include "feCompressedRowStorage.h"

extern int FE_VERBOSE;

bool slepcWasInitialized = false;

void slepcInitialize(int argc, char **argv)
{
#if defined(HAVE_SLEPC)
  if(!slepcWasInitialized) {
    slepcWasInitialized = true;
    PetscErrorCode ierr = SlepcInitialize(&argc, &argv, (char *)0, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    printf("In slepcInitialize : Error : SLEPc was already initialized\n");
    return;
  }
#endif
}

void slepcFinalize()
{
#if defined(HAVE_SLEPC)
  SlepcFinalize();
#endif
}

feStatus createEigenProblem(feEigenProblem *&eigenProblem,
                            eigenSolverType type,
                            std::vector<feSpace *> allFESpaces,
                            std::vector<feBilinearForm *> lhsForms,
                            std::vector<feBilinearForm *> rhsForms,
                            feMetaNumber *metaNumber,
                            feMesh *mesh, int argc, char **argv)
{
  if(metaNumber->getNbUnknowns() == 0)
    return feErrorMsg(FE_STATUS_ERROR, "0 unknowns : attempting to create a 0 x 0 linear system.");

  // Check that all bilinear forms are defined on existing connectivities.
  // The finite element spaces of a single form must be defined on the same geometric connectivity.
  for(feBilinearForm *form : lhsForms) {
    bool OK = false;
    for(feSpace *space : allFESpaces) {
      if(space->getCncGeoTag() == form->getCncGeoTag()) OK = true;
    }
    if(!OK) {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "Bilinear form %s is not defined on any known geometric connectivity (named subdomain).",
        form->getWeakFormName().c_str());
    }
  }
  for(feBilinearForm *form : rhsForms) {
    bool OK = false;
    for(feSpace *space : allFESpaces) {
      if(space->getCncGeoTag() == form->getCncGeoTag()) OK = true;
    }
    if(!OK) {
      return feErrorMsg(
        FE_STATUS_ERROR,
        "Bilinear form %s is not defined on any known geometric connectivity (named subdomain).",
        form->getWeakFormName().c_str());
    }
  }

  switch(type) {
    case SLEPC:
#if defined(HAVE_SLEPC)
      if(argc == 0 || argv == nullptr) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Please provide argc and argv to create a SLEPc eigen solver.");
      }
      if(!slepcWasInitialized) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "SLEPc was not initialized : please initialize SLEPc first by calling "
                          "slepcInitialize(argc, argv) as the very first line of your program.");
      }
      eigenProblem = new feEigenProblem(argc, argv, lhsForms, rhsForms, metaNumber, mesh);
      break;
#else
      return feErrorMsg(FE_STATUS_ERROR,
                        "feNG must be compiled with MPI, PETSc and SLEPc to solve with SLEPc.");
#endif
    default:
      return feErrorMsg(FE_STATUS_ERROR, "Unsupported eigen solver type.");
  }
  return FE_STATUS_OK;
}

feEigenProblem::feEigenProblem(int argc, char **argv,
                 std::vector<feBilinearForm *> lhsForms,
                 std::vector<feBilinearForm *> rhsForms,
                 feMetaNumber *metaNumber,
                 feMesh *mesh)
  : _argc(argc), _argv(argv), _metaNumber(metaNumber), _mesh(mesh)
#if defined(HAVE_SLEPC)
  , _nInc(metaNumber->getNbUnknowns()), _nDofs(metaNumber->getNbDOFs())
#endif
 {
  _numAForms = lhsForms.size();
  _numBForms = rhsForms.size();

#if defined(HAVE_OMP)
    int nThreads = omp_get_max_threads();
#else
    int nThreads = 1;
#endif

  // Copy the bilinear forms    
  for(int i = 0; i < nThreads; ++i) {
    for(feBilinearForm *f : lhsForms) {
      #if defined(HAVE_OMP)
          feBilinearForm *fCpy = new feBilinearForm(*f);
          _AForms.push_back(fCpy);
      #else
          _AForms.push_back(f);
      #endif
    }
    for(feBilinearForm *f : rhsForms) {
      #if defined(HAVE_OMP)
          feBilinearForm *fCpy = new feBilinearForm(*f);
          _BForms.push_back(fCpy);
      #else
          _BForms.push_back(f);
      #endif
    }
  }

 #if defined(HAVE_SLEPC)
  PetscErrorCode ierr;

  // Allocate matrices:
  // Determine the nonzero structure for A and B
  feCompressedRowStorage CRS_A(_metaNumber, _mesh, _AForms, _numAForms);
  feCompressedRowStorage CRS_B(_metaNumber, _mesh, _BForms, _numBForms);
  feInt *NNZ_A = CRS_A.getNnz();
  feInt *NNZ_B = CRS_B.getNnz();
  std::vector<PetscInt> nnz_A(_nInc, 0);
  std::vector<PetscInt> nnz_B(_nInc, 0);
  for(int i = 0; i < _nInc; ++i) {
    nnz_A[i] = NNZ_A[i];
    nnz_B[i] = NNZ_B[i];
  }

  ierr = MatCreate(PETSC_COMM_WORLD, &_A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatCreate(PETSC_COMM_WORLD, &_B); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetSizes(_B, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(_A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(_B); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSeqAIJSetPreallocation(_A, 0, nnz_A.data()); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSeqAIJSetPreallocation(_B, 0, nnz_B.data()); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetUp(_A);
  ierr = MatSetUp(_B);

  PetscInt Ma, Na, Mb, Nb;
  MatGetSize(_A, &Ma, &Na);
  MatGetSize(_B, &Mb, &Nb);
  feInfoCond(FE_VERBOSE > 0, "Created an eigen problem with matrix A of size %d x %d", Ma, Na);
  feInfoCond(FE_VERBOSE > 0, "                         with matrix B of size %d x %d", Mb, Nb);
 #endif
};

void feEigenProblem::viewLHSMatrix()
{
#if defined(HAVE_SLEPC)
  MatView(_A, PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void feEigenProblem::viewRHSMatrix()
{
#if defined(HAVE_SLEPC)
  MatView(_B, PETSC_VIEWER_STDOUT_WORLD);
#endif
}

void feEigenProblem::setMatricesToZero()
{
#if defined(HAVE_SLEPC)
  PetscErrorCode ierr;
  ierr = MatZeroEntries(_A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatZeroEntries(_B); CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feEigenProblem::assembleLHSMatrix(feSolution *sol)
{
#if defined(HAVE_SLEPC)
  PetscErrorCode ierr = 0;

  for(feInt eq = 0; eq < _numAForms; ++eq) {
    feBilinearForm *f = _AForms[eq];
    feCncGeo *cnc = f->getCncGeo();
    int nbColor = cnc->getNbColor();
    std::vector<int> &nbElmPerColor = cnc->getNbElmPerColor();
    std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int nbElmC;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < nbColor; ++iColor) {
      nbElmC = nbElmPerColor[iColor]; // nbElm : nombre d'elm de meme couleur
      listElmC = listElmPerColor[iColor];

      int eqt, elm, numThread = 0;
      double **Ae;
      feInt sizeI, sizeJ;
      std::vector<feInt> niElm;
      std::vector<feInt> njElm;
      std::vector<feInt> adrI;
      std::vector<feInt> adrJ;
      std::vector<PetscScalar> values;

#if defined(HAVE_OMP)
#pragma omp parallel for private(numThread, elm, eqt, f, niElm, njElm, sizeI, sizeJ, adrI, adrJ,   \
                               Ae, ierr, values)
#endif
      for(int iElm = 0; iElm < nbElmC; ++iElm) {
#if defined(HAVE_OMP)
        numThread = omp_get_thread_num();
        eqt = eq + numThread * _numAForms;
        f = _AForms[eqt];
#endif
        elm = listElmC[iElm];
        f->computeMatrix(_metaNumber, _mesh, sol, elm);
        Ae = f->getAe();

        // Determine assignment indices
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        sizeI = adrI.size();
        niElm.reserve(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          if(adrI[i] < _nInc) niElm.push_back(i);
        }
        sizeJ = adrJ.size();
        njElm.reserve(sizeJ);
        for(feInt i = 0; i < sizeJ; ++i) {
          if(adrJ[i] < _nInc) njElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrI.end());
        adrJ.erase(std::remove_if(adrJ.begin(), adrJ.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrJ.end());

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();
        values.resize(sizeI * sizeJ);
        for(feInt i = 0; i < sizeI; ++i) {
          for(feInt j = 0; j < sizeJ; ++j) {
            values[sizeI * i + j] = Ae[niElm[i]][njElm[j]];
          }
        }
        ierr = MatSetValues(_A, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(),
                            ADD_VALUES);
        niElm.clear();
        njElm.clear();
      }
    } // nbColor
  } // NumberOfBilinearForms

  ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Write matrix to a Matlab file
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_A.m", &viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(_A,viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);
#endif
}

void feEigenProblem::assembleRHSMatrix(feSolution *sol)
{
#if defined(HAVE_SLEPC)
  PetscErrorCode ierr = 0;

  for(feInt eq = 0; eq < _numBForms; ++eq) {
    feBilinearForm *f = _BForms[eq];
    feCncGeo *cnc = f->getCncGeo();
    int nbColor = cnc->getNbColor();
    std::vector<int> &nbElmPerColor = cnc->getNbElmPerColor();
    std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int nbElmC;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < nbColor; ++iColor) {
      nbElmC = nbElmPerColor[iColor]; // nbElm : nombre d'elm de meme couleur
      listElmC = listElmPerColor[iColor];

      int eqt, elm, numThread = 0;
      double **Ae;
      feInt sizeI, sizeJ;
      std::vector<feInt> niElm;
      std::vector<feInt> njElm;
      std::vector<feInt> adrI;
      std::vector<feInt> adrJ;
      std::vector<PetscScalar> values;

#if defined(HAVE_OMP)
#pragma omp parallel for private(numThread, elm, eqt, f, niElm, njElm, sizeI, sizeJ, adrI, adrJ,   \
                               Ae, ierr, values)
#endif
      for(int iElm = 0; iElm < nbElmC; ++iElm) {
#if defined(HAVE_OMP)
        numThread = omp_get_thread_num();
        eqt = eq + numThread * _numBForms;
        f = _BForms[eqt];
#endif
        elm = listElmC[iElm];
        f->computeMatrix(_metaNumber, _mesh, sol, elm);
        Ae = f->getAe();

        // Determine assignment indices
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        sizeI = adrI.size();
        niElm.reserve(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          if(adrI[i] < _nInc) niElm.push_back(i);
        }
        sizeJ = adrJ.size();
        njElm.reserve(sizeJ);
        for(feInt i = 0; i < sizeJ; ++i) {
          if(adrJ[i] < _nInc) njElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrI.end());
        adrJ.erase(std::remove_if(adrJ.begin(), adrJ.end(),
                                  [this](const int &x) { return x >= this->_nInc; }),
                   adrJ.end());

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();
        values.resize(sizeI * sizeJ);
        for(feInt i = 0; i < sizeI; ++i) {
          for(feInt j = 0; j < sizeJ; ++j) {
            values[sizeI * i + j] = Ae[niElm[i]][njElm[j]];
          }
        }
        ierr = MatSetValues(_B, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(),
                            ADD_VALUES);
        niElm.clear();
        njElm.clear();
      }
    } // nbColor
  } // NumberOfBilinearForms

  ierr = MatAssemblyBegin(_B, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_B, MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Write matrix to a Matlab file
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_B.m", &viewer);
  PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(_B,viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);
#endif
}

void feEigenProblem::assemble(feSolution *sol)
{
  this->assembleLHSMatrix(sol);
  this->assembleRHSMatrix(sol);
}

void feEigenProblem::solve(feSolution *sol)
{
#if defined (HAVE_SLEPC)

	sol->initializeUnknowns(_mesh, _metaNumber);
  sol->initializeEssentialBC(_mesh, _metaNumber);

	// Create eigensolver context
  EPSCreate(PETSC_COMM_WORLD, &_eps);
  // Set operators. In this case, it is a generalized eigenvalue problem.
  EPSSetOperators(_eps, _A, _B);
  // Solve for the N smallest eigenpairs
  EPSSetWhichEigenpairs(_eps, EPS_SMALLEST_REAL);
  // Set solver parameters at runtime
  EPSSetFromOptions(_eps);

  feInfo("=====================================================================");
  feInfo("                     Information from SLEPC                          ");
  feInfo("=====================================================================");

	PetscInt nev, maxit, i, its, lits, nconv;
	EPSSolve(_eps);
	EPSGetIterationNumber(_eps,&its);
  PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its);
  EPSGetST(_eps,&_st);
  STGetKSP(_st,&_ksp);
  KSPGetTotalIterations(_ksp,&lits);
  PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %" PetscInt_FMT "\n",lits);
  EPSGetType(_eps,&_type);
  PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",_type);
  EPSGetDimensions(_eps,&nev,NULL,NULL);
  PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev);
  EPSGetTolerances(_eps,&_tol,&maxit);
  PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)_tol,maxit);

  EPSErrorView(_eps,EPS_ERROR_RELATIVE,NULL);
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);
  EPSConvergedReasonView(_eps,PETSC_VIEWER_STDOUT_WORLD);
  EPSErrorView(_eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);

  EPSGetConverged(_eps, &nconv);

  feInfo("=====================================================================");
  feInfo("                     End information from SLEPC                      ");
  feInfo("=====================================================================");

  PetscScalar eigr, eigi;

  // Store eigenpairs in the eigenProblem.
  _eigenPairs.clear();
  for (i = 0; i < nconv; i++){
  	_eigenPairs.push_back(eigenPair());
  	MatCreateVecs(_A, NULL, &_eigenPairs[i].vecReal);
  	MatCreateVecs(_A, NULL, &_eigenPairs[i].vecImag);
  	EPSGetEigenpair(_eps, i, &eigr, &eigi, _eigenPairs[i].vecReal, _eigenPairs[i].vecImag);
  	_eigenPairs[i].valReal = eigr;
  	_eigenPairs[i].valImag = eigi;
  }

  // Start with the eigenpair associated to the lowest eigenvalue
  // std::reverse(_eigenPairs.begin(),_eigenPairs.end());

  // Re-normalize by setting the max nodal value to +-1
  PetscScalar norm;
  for (i = 0; i < nconv; i++){
  	VecNorm(_eigenPairs[i].vecReal, NORM_INFINITY, &norm);
  	VecScale(_eigenPairs[i].vecReal, 1./norm);
  }

#endif
}

// Sets the eigenPairIndex-th eigenmode as the active solution in the feSolution.
// For example, if U is the computed field, the eigenmode will be stored in the feSolution instead of U.
void feEigenProblem::setEigenmodeAsActiveSolution(feSolution *sol, size_t eigenPairIndex, std::string fieldID)
{
#if defined(HAVE_SLEPC)
	std::vector<int> allUnknownDOF = _metaNumber->getUnknownDOF(fieldID);

	PetscScalar *vecRealArray;
  VecGetArray(_eigenPairs[eigenPairIndex].vecReal, &vecRealArray);
	for(int dof : allUnknownDOF){
		sol->setSolAtDOF(dof, vecRealArray[dof]);
	}
	VecRestoreArray(_eigenPairs[eigenPairIndex].vecReal, &vecRealArray);
#endif
}

feEigenProblem::~feEigenProblem()
{
#if defined(HAVE_SLEPC)
	for(size_t i = 0; i < _eigenPairs.size(); ++i){
		VecDestroy(&_eigenPairs[i].vecReal);
		VecDestroy(&_eigenPairs[i].vecImag);
	}
#endif
}