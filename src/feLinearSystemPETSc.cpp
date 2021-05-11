#include "feLinearSystemPETSc.h"

void feLinearSystemPETSc::initialize(){
  PetscErrorCode ierr;
  PetscMPIInt    size;

  ierr = PetscInitialize(&_argc, &_argv,(char*) 0, nullptr);
  if(ierr){
    printf("In feLinearSystemPETSc::feLinearSystemPETSc : Error at PetscInitialize()\n");
    return;
  }
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRABORT(PETSC_COMM_WORLD, ierr); //CHKERRMPI(ierr);
  if(size != 1){
    // SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
    printf("In feLinearSystemPETSc::feLinearSystemPETSc : This is a uniprocessor example only\n");
    return;
  }
  // Allocate vectors for the residual _res and the solution correction _dx
  ierr = PetscOptionsGetInt(NULL, NULL, "-n", &_nInc, NULL);    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &_dx);                     CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = PetscObjectSetName((PetscObject) _dx, "Solution");     CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetSizes(_dx, PETSC_DECIDE, _nInc);                 CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSetFromOptions(_dx);                                CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDuplicate(_dx, &_res);                              CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // Allocate the matrix _A
  ierr = MatCreate(PETSC_COMM_WORLD, &_A);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetSizes(_A,PETSC_DECIDE,PETSC_DECIDE,_nInc,_nInc); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(_A);                                 CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetUp(_A);                                          CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // Create the Krylov solver (default is GMRES)
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetOperators(ksp,_A,_A);                            CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // Set preconditioner
  ierr = KSPGetPC(ksp,&preconditioner);                         CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PCSetType(preconditioner,PCJACOBI);                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = PCSetType(preconditioner,PCILU);                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetTolerances(ksp,1.e-16,1e-16,PETSC_DEFAULT,100);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);                                CHKERRABORT(PETSC_COMM_WORLD, ierr);
};

void feLinearSystemPETSc::viewMatrix(){
  MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
  PetscViewer viewer;
  PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,600,600,&viewer);
  PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
  PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
  MatView(_A,viewer);
}

void feLinearSystemPETSc::setToZero(){
  if(recomputeMatrix){
    this->setMatrixToZero();
  }
  this->setResidualToZero();
}

void feLinearSystemPETSc::setMatrixToZero(){
  PetscErrorCode ierr = 0;
  ierr = MatZeroEntries(_A);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // viewMatrix();
}

void feLinearSystemPETSc::setResidualToZero(){
  // viewMatrix();
  VecZeroEntries(_res);
}

void feLinearSystemPETSc::assembleMatrices(feSolution *sol){
  PetscErrorCode ierr = 0;
  if(recomputeMatrix){
    PetscInt I,J;
    int niElm, njElm;
    std::vector<int> adrI, adrJ;

    for(feBilinearForm *f : _formMatrices){
      int nElm = _mesh->getNbElm(f->getCncGeoTag());
      for(int iElm = 0; iElm < nElm; ++iElm){
        f->computeMatrix(_metaNumber, _mesh, sol, iElm); // Matrice elementaire
        std::vector<double> Ae = f->getAe();
        niElm = f->getNiElm();
        njElm = f->getNjElm();
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        for(int i = 0; i < niElm; ++i){
          I = adrI[i];
          if(I < _nInc){
            for(int j = 0; j < njElm; ++j){
              J = adrJ[j];
              if(J < _nInc){
                ierr = MatSetValue(_A, I, J, Ae[njElm*i+j], ADD_VALUES); // TODO : assigner par blocs
                CHKERRABORT(PETSC_COMM_WORLD, ierr);
              }
            }
          }
        }
      }
    }
    ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
  } // if(recomputeMatrix)
  // viewMatrix();
}

void feLinearSystemPETSc::assembleResiduals(feSolution *sol){
  PetscErrorCode ierr;
  PetscInt I;
  int niElm;
  std::vector<int> adrI;
  for(feBilinearForm *f : _formResiduals){
    int nElm = _mesh->getNbElm(f->getCncGeoTag());
    for(int iElm = 0; iElm < nElm; ++iElm){
      f->computeResidual(_metaNumber, _mesh, sol, iElm); // Residu elementaire
      std::vector<double> Be = f->getBe();
      niElm = f->getNiElm();
      adrI = f->getAdrI();
      for(int i = 0; i < niElm; ++i){
        I = adrI[i];
        if(I < _nInc){
          ierr = VecSetValue(_res, I, Be[i], ADD_VALUES); // TODO : assigner par blocs
          CHKERRABORT(PETSC_COMM_WORLD, ierr);
        }
      }
    }
  }
  // VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
}

void feLinearSystemPETSc::assemble(feSolution *sol){
  this->assembleMatrices(sol);
  this->assembleResiduals(sol);
}

void feLinearSystemPETSc::solve(double *normDx, double *normResidual){
  // PetscInt its;
  PetscErrorCode ierr = KSPSolve(ksp,_res,_dx); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  // KSPGetIterationNumber(ksp,&its);
  // PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
  ierr = VecNorm(_res, NORM_2, normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm( _dx, NORM_2, normDx);       CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void feLinearSystemPETSc::correctSolution(feSolution *sol){
  PetscScalar *array;
  VecGetArray(_dx, &array);
  for(int i = 0; i < _nInc; ++i)
    sol->incrementSolAtDOF(i, array[i]);
  VecRestoreArray(_dx, &array);
}

void feLinearSystemPETSc::assignResidualToDCResidual(feSolutionContainer *solContainer){
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i)
    solContainer->_fResidual[0][i] = array[i];
  VecRestoreArray(_res, &array);
}

void feLinearSystemPETSc::printResidual(){
  // VecView(_res, PETSC_VIEWER_STDOUT_WORLD);
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i)
    printf("%12.16f\n", array[i]);
  VecRestoreArray(_res, &array);
}

// void feLinearSystemPETSc::check(){
//   VecView(_dx,PETSC_VIEWER_STDOUT_WORLD);

//   PetscReal normRes, normDx;
//   PetscErrorCode ierr = VecNorm(_res,NORM_2,&normRes); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   ierr = VecNorm( _dx,NORM_2,&normDx); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   std::cout<<"Norme res = "<<normRes<<std::endl;
//   std::cout<<"Norme dx = "<<normDx<<std::endl;
//   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   // ierr = VecAXPY(x,-1.0,u);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   // ierr = VecNorm(x,NORM_2,&norm);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   // ierr = KSPGetIterationNumber(ksp,&its);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRABORT(PETSC_COMM_WORLD, ierr);
// }
  
void feLinearSystemPETSc::finalize(){
  PetscErrorCode ierr;
  ierr = KSPDestroy(&ksp);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatDestroy(&_A);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecDestroy(&_res); CHKERRABORT(PETSC_COMM_WORLD, ierr); 
  ierr = VecDestroy(&_dx);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PetscFinalize();
}