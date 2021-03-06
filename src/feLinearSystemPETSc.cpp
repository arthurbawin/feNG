#include "feLinearSystemPETSc.h"
#include "feCompressedRowStorage.h"

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
    printf("In petscInitialize : Error : PETSc was already initialized\n");
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
  PetscErrorCode ierr;
  PetscMPIInt size;

  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  if(size != 1) {
    printf("In feLinearSystemPETSc::initialize : This is a uniprocessor example only\n");
    return;
  }
  // Allocate vectors and matrix
  ierr = VecCreate(PETSC_COMM_WORLD, &_dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetSizes(_dx, PETSC_DECIDE, _nInc);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(_dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = VecDuplicate(_dx, &_res);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(_dx, &_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSet(_dx, 1.0);

  // Determine the nonzero structure
  feCompressedRowStorage CRS(_metaNumber, _mesh, _formMatrices);
  feInt *NNZ = CRS.getNnz();
  std::vector<PetscInt> nnz(_nInc, 0);
  for(int i = 0; i < _nInc; ++i) {
    nnz[i] = NNZ[i];
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

  PetscReal rel_tol = 1e-6;
  PetscReal abs_tol = 1e-12;
  PetscReal div_tol = 1e6;
  PetscInt max_iter = 500;

  ierr = KSPSetTolerances(ksp, rel_tol, abs_tol, div_tol, max_iter);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  PetscInt M, N, Nres, Ndx, NlinSysRes;
  MatGetSize(_A, &M, &N);
  VecGetSize(_res, &Nres);
  VecGetSize(_dx, &Ndx);
  VecGetSize(_linSysRes, &NlinSysRes);
  feInfoCond(FE_VERBOSE > 0, "Created a linear system of size %d x %d\n", M, N);
  feInfoCond(FE_VERBOSE > 0, "Created a res vector of size %d\n", Nres);
  feInfoCond(FE_VERBOSE > 0, "Created a dx vector of size %d\n", Ndx);
  feInfoCond(FE_VERBOSE > 0, "Created a linSysRes vector of size %d\n", NlinSysRes);
#endif
}

feLinearSystemPETSc::feLinearSystemPETSc(int argc, char **argv,
                                         std::vector<feBilinearForm *> bilinearForms,
                                         feMetaNumber *metaNumber, feMesh *mesh)
  : feLinearSystem(bilinearForms, metaNumber, mesh), _argc(argc), _argv(argv)
#if defined(HAVE_PETSC)
    ,
    _nInc(metaNumber->getNbUnknowns()), _nDofs(metaNumber->getNbDOFs())
#endif
{
  this->initialize();
}

void feLinearSystemPETSc::viewMatrix()
{
#if defined(HAVE_PETSC)
  MatView(_A, PETSC_VIEWER_STDOUT_WORLD);
  // PetscViewer viewer;
  // PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, NULL, 0, 0, 600, 600, &viewer);
  // PetscObjectSetName((PetscObject)viewer, "Line graph Plot");
  // PetscViewerPushFormat(viewer, PETSC_VIEWER_DRAW_LG);
  // MatView(_A, viewer);
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
  PetscErrorCode ierr = 0;
  ierr = MatZeroEntries(_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // viewMatrix();
#endif
}

void feLinearSystemPETSc::setResidualToZero()
{
#if defined(HAVE_PETSC)
  VecZeroEntries(_res);
#endif
}

void feLinearSystemPETSc::assembleMatrices(feSolution *sol)
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = 0;
  if(recomputeMatrix) {
    // tic();
    PetscInt I, J;
    std::vector<PetscScalar> values;
    int sizeI, sizeJ;
    std::vector<int> niElm;
    std::vector<int> njElm;
    std::vector<int> adrI;
    std::vector<int> adrJ;

    for(feBilinearForm *f : _formMatrices) {
      int nElm = _mesh->getNbElm(f->getCncGeoTag());
      for(int iElm = 0; iElm < nElm; ++iElm) {
        f->computeMatrix(_metaNumber, _mesh, sol, iElm); // Matrice elementaire
        double **Ae = f->getAe();
        // Determine assignment indices
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        sizeI = adrI.size();
        niElm.reserve(sizeI);
        for(int i = 0; i < sizeI; ++i) {
          if(adrI[i] < _nInc) niElm.push_back(i);
        }
        sizeJ = adrJ.size();
        njElm.reserve(sizeJ);
        for(int i = 0; i < sizeJ; ++i) {
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
        for(int i = 0; i < sizeI; ++i) {
          for(int j = 0; j < sizeJ; ++j) {
            values[sizeI * i + j] = Ae[niElm[i]][njElm[j]];
          }
        }
        ierr = MatSetValues(_A, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(),
                            ADD_VALUES);
        niElm.clear();
        njElm.clear();
      }
    }
    ierr = MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // double normMat = 0.0;
    // ierr = MatNorm(_A, NORM_FROBENIUS, &normMat); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // printf("Norme de la matrice : %10.10e\n", normMat);
    // toc();
  } // if(recomputeMatrix)
#endif
}

void feLinearSystemPETSc::assembleResiduals(feSolution *sol)
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr;
  std::vector<PetscScalar> values;
  int sizeI;
  std::vector<int> niElm;
  std::vector<int> adrI;
  double normResidual = 0.0;
  ierr = VecNorm(_res, NORM_2, &normResidual);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  for(feBilinearForm *f : _formResiduals) {
    int nElm = _mesh->getNbElm(f->getCncGeoTag());
    for(int iElm = 0; iElm < nElm; ++iElm) {
      f->computeResidual(_metaNumber, _mesh, sol, iElm); // Residu elementaire
      double *Be = f->getBe();
      // Determine assignment indices
      adrI = f->getAdrI();
      sizeI = adrI.size();
      niElm.reserve(sizeI);
      for(int i = 0; i < sizeI; ++i) {
        if(adrI[i] < _nInc) niElm.push_back(i);
      }
      adrI.erase(
        std::remove_if(adrI.begin(), adrI.end(), [this](const int &x) { return x >= this->_nInc; }),
        adrI.end());
      // Copy Be into vector
      sizeI = adrI.size();
      values.resize(sizeI);
      for(int i = 0; i < sizeI; ++i) {
        values[i] = Be[niElm[i]];
        // std::cout<< "les valeurs dans Be sont "<< values[i] << std::endl;
      }
      ierr = VecSetValues(_res, adrI.size(), adrI.data(), values.data(), ADD_VALUES);
      niElm.clear();
    }
  }
  // VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
  // normResidual = 0.0;
  // ierr = VecNorm(_res, NORM_2, &normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // printf("Norme du r??sidu : %10.10e\n", normResidual);
#endif
}

void feLinearSystemPETSc::assemble(feSolution *sol)
{
  this->assembleMatrices(sol);
  this->assembleResiduals(sol);
}

// Solve the system and compute norms of solution and residuals
void feLinearSystemPETSc::solve(double *normDx, double *normResidual, double *normAxb, int *nIter)
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr = KSPSolve(ksp, _res, _dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  KSPGetIterationNumber(ksp, nIter);
  VecSet(_linSysRes, 0.0);
  MatMult(_A, _dx, _linSysRes);
  VecAXPY(_linSysRes, -1.0, _res);
  ierr = VecNorm(_res, NORM_MAX, normResidual);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_dx, NORM_MAX, normDx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm(_linSysRes, NORM_MAX, normAxb);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
#endif
}

void feLinearSystemPETSc::correctSolution(feSolution *sol)
{
#if defined(HAVE_PETSC)
  std::vector<double> &_sol = sol->getSolutionReference();
  PetscScalar *array;
  VecGetArray(_dx, &array);
  for(int i = 0; i < _nInc; ++i) _sol[i] += array[i];
  VecRestoreArray(_dx, &array);
#endif
}

void feLinearSystemPETSc::assignResidualToDCResidual(feSolutionContainer *solContainer)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i) solContainer->_fResidual[0][i] = array[i];
  VecRestoreArray(_res, &array);
#endif
}

void feLinearSystemPETSc::applyCorrectionToResidual(double coeff, std::vector<double> &d)
{
#if defined(HAVE_PETSC)
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i) array[i] += coeff * d[i];
  VecRestoreArray(_res, &array);
#endif
}

void feLinearSystemPETSc::printResidual()
{
#if defined(HAVE_PETSC)
  // VecView(_res, PETSC_VIEWER_STDOUT_WORLD);
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i) printf("%12.16f\n", array[i]);
  VecRestoreArray(_res, &array);
#endif
}

void feLinearSystemPETSc::finalize()
{
#if defined(HAVE_PETSC)
  PetscErrorCode ierr;
  ierr = KSPDestroy(&ksp);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  printf("In feLinearSystemPETSc::feLinearSystemPETSc : Deleted KSP\n");
  ierr = MatDestroy(&_A);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  printf("In feLinearSystemPETSc::feLinearSystemPETSc : Deleted matrix A         of size %d\n",
         _nInc);
  ierr = VecDestroy(&_linSysRes);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  printf("In feLinearSystemPETSc::feLinearSystemPETSc : Deleted vector linSysRes of size %d\n",
         _nInc);
  ierr = VecDestroy(&_res);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  printf("In feLinearSystemPETSc::feLinearSystemPETSc : Deleted vector res       of size %d\n",
         _nInc);
  ierr = VecDestroy(&_dx);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  printf("In feLinearSystemPETSc::feLinearSystemPETSc : Deleted vector dx        of size %d\n",
         _nInc);
#endif
}

feLinearSystemPETSc::~feLinearSystemPETSc() { this->finalize(); }
