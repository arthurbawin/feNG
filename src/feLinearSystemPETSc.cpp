#include "feLinearSystemPETSc.h"
#include "feCompressedRowStorage.h"

#if defined(HAVE_PETSC)
#include "petscdraw.h"
#endif

extern int FE_VERBOSE;

bool petscWasInitialized = false;

void petscInitialize(int argc, char **argv)
{
#if defined(HAVE_PETSC)
  if(!petscWasInitialized) {
    petscWasInitialized = true;
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, (char *)0, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } else {
    feErrorMsg(FE_STATUS_ERROR, "PETSc was already initialized\n");
    return;
  }
#endif
}

void petscFinalize()
{
#if defined(HAVE_PETSC)
  PetscFinalize();
#endif
}

void feLinearSystemPETSc::initialize()
{
#if defined(HAVE_PETSC)
  PetscMPIInt size;
  PetscErrorCode ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  if(size != 1) {
    feErrorMsg(FE_STATUS_ERROR, "Cannot run on more than 1 processor :/");
    exit(-1);
  }

  // Allocate vectors and matrix
  ierr = VecCreate(PETSC_COMM_WORLD, &_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetSizes(_du, PETSC_DECIDE, _nInc);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecDuplicate(_du, &_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(_du, &_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSet(_du, 1.0);

  // Determine the nonzero structure
  // feCompressedRowStorage CRS(_nInc, _formMatrices, _numMatrixForms);
  _EZCRS = new feEZCompressedRowStorage(_nInc, _formMatrices, _numMatrixForms);

  // feInt *NNZ = CRS.getNnz();
  std::vector<PetscInt> nnz(_nInc, 0);
  for(int i = 0; i < _nInc; ++i) {
    // nnz[i] = NNZ[i];
    nnz[i] = _EZCRS->getNnzAt(i);
  }

  bool withPrealloc = true;
  bool printInfos = false;

  if(withPrealloc) {
    ierr = MatCreate(PETSC_COMM_WORLD, &_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetFromOptions(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSeqAIJSetPreallocation(_A, 0, nnz.data());
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetUp(_A);

    if(printInfos) {
      MatInfo info;
      double mal, nz_a, nz_u, mem, nz_un;

      MatGetInfo(_A, MAT_LOCAL, &info);
      mal = info.mallocs;
      nz_a = info.nz_allocated;
      nz_u = info.nz_used;
      nz_un = info.nz_unneeded;
      mem = info.memory;

      std::cout << "mal = " << mal << std::endl;
      std::cout << "mem = " << mem << std::endl;
      std::cout << "nz_a = " << nz_a << std::endl;
      std::cout << "nz_u = " << nz_u << std::endl;
      std::cout << "nz_un = " << nz_un << std::endl;
    }
  } else {
    // Without allocation (bad) :
    ierr = MatCreate(PETSC_COMM_WORLD, &_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetFromOptions(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatSetUp(_A);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  }

  // Create the Krylov solver (default is GMRES)
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetType(ksp, KSPBCGS);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetOperators(ksp, _A, _A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // Set preconditioner
  ierr = KSPGetPC(ksp, &preconditioner);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = PCSetType(preconditioner, PCILU);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = KSPSetTolerances(ksp, (PetscReal) _rel_tol, (PetscReal) _abs_tol, 
    (PetscReal) _div_tol, (PetscInt) _max_iter);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  PetscInt M, N, Nres, Ndx, NlinSysRes;
  MatGetSize(_A, &M, &N);
  VecGetSize(_rhs, &Nres);
  VecGetSize(_du, &Ndx);
  VecGetSize(_linSysRes, &NlinSysRes);

  feInfoCond(FE_VERBOSE > 0, "\t\tCreated a PETSc linear system of size %d x %d", M, N);
  feInfoCond(FE_VERBOSE > 0, "\t\tIterative linear solver info:");

  KSPType krylovType;
  ierr = KSPGetType(ksp, &krylovType);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tKrylov method: %s", krylovType);

  PCType preconditionerType;
  ierr = PCGetType(preconditioner, &preconditionerType);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tPreconditioner: %s", preconditionerType);

  feInfoCond(FE_VERBOSE > 0, "\t\t\tRelative tolerance: %1.4e", _rel_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tAbsolute tolerance: %1.4e", _abs_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tDivergence criterion: %1.4e", _div_tol);
  feInfoCond(FE_VERBOSE > 0, "\t\t\tMax number of iterations: %d", _max_iter);
#endif
}

feLinearSystemPETSc::feLinearSystemPETSc(int argc, char **argv,
                                         std::vector<feBilinearForm *> bilinearForms,
                                         int numUnknowns)
  : feLinearSystem(bilinearForms), _argc(argc), _argv(argv)
#if defined(HAVE_PETSC)
    , _nInc(numUnknowns)
#endif
{
  this->initialize();
}

void feLinearSystemPETSc::viewMatrix()
{
#if defined(HAVE_PETSC)
  if(_displayMatrixInConsole){
    MatView(_A, PETSC_VIEWER_STDOUT_WORLD);
  }

  if(_displayMatrixInWindow){

    PetscErrorCode ierr = 0;
    // ierr = MatGetOrdering(_A, MATORDERINGRCM, &_rowMap, &_colMap); CHKERRV(ierr);
    // ierr = MatPermute(_A, _rowMap, _colMap, &_Ap); CHKERRV(ierr);

    PetscViewer viewer;
    PetscDraw draw;
    ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 600, 600, &viewer); CHKERRV(ierr);
    ierr = MatView(_A, viewer); CHKERRV(ierr);
    ierr = PetscViewerDrawGetDraw(viewer, 0, &draw); CHKERRV(ierr);
    ierr = PetscDrawSetPause(draw, -1); CHKERRV(ierr); // Wait for user
    PetscDrawPause(draw);
    ierr = PetscViewerDestroy(&viewer); CHKERRV(ierr);
  }

#endif
}

void feLinearSystemPETSc::viewRHS()
{
#if defined(HAVE_PETSC)
  if(_displayRHSInConsole){
    VecView(_rhs, PETSC_VIEWER_STDOUT_WORLD);
  }
#endif
}

void feLinearSystemPETSc::setToZero()
{
  if(recomputeMatrix) {
    this->setMatrixToZero();
  }
  this->setResidualToZero();
}

void feLinearSystemPETSc::setMatrixToZero()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = MatZeroEntries(_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::setResidualToZero()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = VecZeroEntries(_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::assembleMatrices(feSolution *sol)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  for(feInt eq = 0; eq < _numMatrixForms; ++eq) {

    feBilinearForm *f = _formMatrices[eq];
    feCncGeo *cnc = f->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &numElemPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int numElementsInColor;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor) {
      numElementsInColor = numElemPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      int elm;
      double **Ae;
      feInt sizeI;
      feInt sizeJ;
      std::vector<feInt> niElm;
      std::vector<feInt> njElm;
      std::vector<feInt> adrI;
      std::vector<feInt> adrJ;
      std::vector<PetscScalar> values;

      #if defined(HAVE_OMP)
      #pragma omp parallel for private(elm, f, niElm, njElm, sizeI, sizeJ, adrI, adrJ,   \
                               Ae, ierr, values)
      #endif
      for(int iElm = 0; iElm < numElementsInColor; ++iElm) {
        #if defined(HAVE_OMP)
        int eqt = eq + omp_get_thread_num() * _numMatrixForms;
        f = _formMatrices[eqt];
        #endif
        elm = listElmC[iElm];

        // Compute element-wise matrix
        f->computeMatrix(sol, elm);
        Ae = f->getAe();

        // Determine global assignment indices
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
          [this](const int &x) { return x >= this->_nInc; }), adrI.end());
        adrJ.erase(std::remove_if(adrJ.begin(), adrJ.end(),
          [this](const int &x) { return x >= this->_nInc; }), adrJ.end());

        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();

        values.resize(sizeI * sizeJ);
        for(feInt i = 0; i < sizeI; ++i) {
          for(feInt j = 0; j < sizeJ; ++j) {
            values[sizeJ * i + j] = Ae[niElm[i]][njElm[j]];
          }
        }

        // Increment global matrix
        ierr = MatSetValues(_A, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(), ADD_VALUES);
        niElm.clear();
        njElm.clear();
      }
    }
  }

  ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // double normMat;
  // ierr = MatNorm(_A, NORM_FROBENIUS, &normMat);
  // CHKERRABORT(PETSC_COMM_WORLD, ierr);

  viewMatrix();

  if(_exportMatrixMatlab){
    feInfoCond(FE_VERBOSE > 1, "\t\t\tExporting global matrix to Matlab");
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix.m", &viewer);
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    MatView(_A,viewer);
    PetscViewerPopFormat(viewer);
    PetscViewerDestroy(&viewer);
  }
  
#endif
}

void feLinearSystemPETSc::assembleResiduals(feSolution *sol)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  for(feInt eq = 0; eq < _numResidualForms; ++eq) {

    feBilinearForm *f = _formResiduals[eq];
    feCncGeo *cnc = f->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &numElemPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int numElementsInColor;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor) {
      int numElementsInColor = numElemPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      int elm;
      double *Be;
      feInt sizeI;
      std::vector<feInt> niElm;
      std::vector<feInt> adrI;
      std::vector<PetscScalar> values;

      #if defined(HAVE_OMP)
      #pragma omp parallel for private(elm, niElm, Be, f, adrI, ierr, values, sizeI)
      #endif
      for(int iElm = 0; iElm < numElementsInColor; ++iElm) {
        #if defined(HAVE_OMP)
        int eqt = eq + omp_get_thread_num() * _numResidualForms;
        f = _formResiduals[eqt];
        #endif

        elm = listElmC[iElm];

        // Compute the element-wise residual
        f->computeResidual(sol, elm);
        Be = f->getBe();

        // Determine global assignment indices
        adrI = f->getAdrI();
        sizeI = adrI.size();
        niElm.reserve(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          if(adrI[i] < _nInc) niElm.push_back(i);
        }

        adrI.erase(std::remove_if(adrI.begin(), adrI.end(), 
          [this](const int &x) { return x >= this->_nInc; }), adrI.end());

        sizeI = adrI.size();
        values.resize(sizeI);
        for(feInt i = 0; i < sizeI; ++i) {
          values[i] = Be[niElm[i]];
        }

        // Increment global residual
        ierr = VecSetValues(_rhs, adrI.size(), adrI.data(), values.data(), ADD_VALUES);
        niElm.clear();
      }
    }
  }

  ierr = VecAssemblyBegin(_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecAssemblyEnd(_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  
  viewRHS();

#endif
}

void feLinearSystemPETSc::assemble(feSolution *sol)
{
  if(recomputeMatrix){
    this->assembleMatrices(sol);
  }
  this->assembleResiduals(sol);
}

// Solve the linear system and compute max norms of solution and residuals
bool feLinearSystemPETSc::solve(double *normSolution, double *normRHS, double *normResidualAxMinusb, int *nIter)
{
#if defined(HAVE_PETSC)

  PetscErrorCode ierr = 0;

  // Set tolerances, in case the have changed since creation
  ierr = KSPSetTolerances(ksp, (PetscReal) _rel_tol, (PetscReal) _abs_tol, 
    (PetscReal) _div_tol, (PetscInt) _max_iter);

  // Solve
  ierr = KSPSolve(ksp, _rhs, _du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason);

  if(reason < 0){
    // Solve failed
    const char *reasonString;
    KSPGetConvergedReasonString(ksp, &reasonString);
    feWarning("\nPETSc solve failed with KSP error code: %s", reasonString);
    if(reason == KSP_DIVERGED_PC_FAILED){
      PCFailedReason reasonpc;
      ierr = PCGetFailedReason(preconditioner, &reasonpc);
      switch(reasonpc){
        case PC_SETUP_ERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_SETUP_ERROR"); break;
        case PC_NOERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_NOERROR"); break;
        case PC_FACTOR_STRUCT_ZEROPIVOT:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_FACTOR_STRUCT_ZEROPIVOT"); break;
        case PC_FACTOR_NUMERIC_ZEROPIVOT:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_FACTOR_NUMERIC_ZEROPIVOT"); break;
        case PC_FACTOR_OUTMEMORY:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_FACTOR_OUTMEMORY"); break;
        case PC_FACTOR_OTHER:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_FACTOR_OTHER"); break;
        case PC_INCONSISTENT_RHS:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_INCONSISTENT_RHS"); break;
        case PC_SUBPC_ERROR:
          feWarning("\nPETSc solve failed with preconditioner error code: %s", "PC_SUBPC_ERROR"); break;
      }
    }
  }

  KSPGetIterationNumber(ksp, nIter);
  VecSet(_linSysRes, 0.0);
  MatMult(_A, _du, _linSysRes);
  VecAXPY(_linSysRes, -1.0, _rhs);
  ierr = VecNorm(_rhs,       NORM_MAX, normRHS); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_du,        NORM_MAX, normSolution); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_linSysRes, NORM_MAX, normResidualAxMinusb); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  return reason > 0;
#endif
  return false;
}

void feLinearSystemPETSc::correctSolution(feSolution *sol)
{
#if defined(HAVE_PETSC)
  std::vector<double> &_sol = sol->getSolutionReference();
  PetscScalar *array;
  VecGetArray(_du, &array);
  for(int i = 0; i < _nInc; ++i) _sol[i] += array[i];
  VecRestoreArray(_du, &array);
#endif
}

void feLinearSystemPETSc::assignResidualToDCResidual(feSolutionContainer *solContainer)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_rhs, &array);
  for(int i = 0; i < _nInc; ++i) solContainer->_fResidual[0][i] = array[i];
  VecRestoreArray(_rhs, &array);
#endif
}

void feLinearSystemPETSc::applyCorrectionToResidual(double coeff, std::vector<double> &d)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_rhs, &array);
  for(int i = 0; i < _nInc; ++i) array[i] += coeff * d[i];
  VecRestoreArray(_rhs, &array);
#endif
}

void feLinearSystemPETSc::finalize()
{
#if defined(HAVE_PETSC)
  delete _EZCRS;
  PetscErrorCode ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatDestroy(&_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_rhs);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_du);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

feLinearSystemPETSc::~feLinearSystemPETSc() { this->finalize(); }
