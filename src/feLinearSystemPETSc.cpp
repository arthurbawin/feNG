#include "feLinearSystemPETSc.h"
#include "feCompressedRowStorage.h"

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
  ierr = VecDuplicate(_dx, &_foo);                              CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecSet(_dx, 1.0);

  // Determine the nonzero structure
  feCompressedRowStorage CRS(_metaNumber, _mesh, _formMatrices);
  feInt *NNZ = CRS.getNnz();
  std::vector<PetscInt> nnz(_nInc,0);
  int sum = 0;
  for(int i = 0; i < _nInc; ++i){
    nnz[i] = NNZ[i];
    sum += NNZ[i];
  }

  // for(auto val : nnz)
  //   std::cout<<val<<std::endl;

  // Allocate the matrix _A
  // ierr = MatCreate(PETSC_COMM_WORLD, &_A);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,_nInc,_nInc, 0, nnz.data(), &_A);   CHKERRABORT(PETSC_COMM_WORLD, ierr);

  // ierr = MatCreate(PETSC_COMM_WORLD, &_A);                               CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetType(_A,MATSEQAIJ);                                        CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, _nInc, _nInc);     CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetOption(_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetFromOptions(_A);                                         CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetOption(_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSeqAIJSetPreallocation(_A, 0, nnz.data());                  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSeqAIJSetPreallocation(_A, 0, nullptr);                  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetOption(_A, MAT_FORCE_DIAGONAL_ENTRIES , PETSC_TRUE);     CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetUp(_A);

  // MatInfo info;
  // double  mal, nz_a, nz_u, mem, nz_un;

  // MatGetInfo(_A,MAT_LOCAL,&info);
  // mal  = info.mallocs;
  // nz_a = info.nz_allocated;
  // nz_u = info.nz_used;
  // nz_un = info.nz_unneeded;
  // mem = info.memory;

  // std::cout<<"Initial sum = "<<sum<<std::endl;
  // std::cout<<"mal = "<<mal<<std::endl;
  // std::cout<<"mem = "<<mal<<std::endl;
  // std::cout<<"nz_a = "<<nz_a<<std::endl;
  // std::cout<<"nz_u = "<<nz_u<<std::endl;
  // std::cout<<"nz_un = "<<nz_un<<std::endl;

  // ierr = MatCreate(PETSC_COMM_WORLD, &_A);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetSizes(_A,PETSC_DECIDE,PETSC_DECIDE,_nInc,_nInc); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetFromOptions(_A);                                 CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = MatSetUp(_A);                                          CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // Create the Krylov solver (default is GMRES)
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);                      CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetType(ksp, KSPBCGS);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetOperators(ksp,_A,_A);                            CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // Set preconditioner
  ierr = KSPGetPC(ksp,&preconditioner);                         CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PCSetType(preconditioner,PCJACOBI);                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = PCSetType(preconditioner,PCILU);                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // ierr = PCSetType(preconditioner,PCNONE);                    CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetTolerances(ksp,1e-12,1e-12,100000,5000);  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = KSPSetFromOptions(ksp);                                CHKERRABORT(PETSC_COMM_WORLD, ierr);
  double normResidual = 0.0;
  double normDx = 0.0;
  ierr = VecNorm(_res, NORM_MAX, &normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm( _dx, NORM_MAX, &normDx);       CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // std::cout<<"Taille du systeme : "<<_nInc<<std::endl;
  // std::cout<<"Avant résolution : norme du residu :        "<<normResidual<<std::endl;
  // std::cout<<"Avant résolution : norme de la correction : "<<normDx<<std::endl;
}

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

// void feLinearSystemPETSc::assembleMatrices(feSolution *sol){
//   PetscErrorCode ierr = 0;
//   if(recomputeMatrix){
//     PetscInt I,J;
//     int niElm, njElm;
//     std::vector<int> adrI, adrJ;

//     for(feBilinearForm *f : _formMatrices){
//       int nElm = _mesh->getNbElm(f->getCncGeoTag());
//       for(int iElm = 0; iElm < nElm; ++iElm){
//         f->computeMatrix(_metaNumber, _mesh, sol, iElm); // Matrice elementaire
//         // std::vector<double> Ae = f->getAe();
//         double** Ae = f->getAe();
//         niElm = f->getNiElm();
//         njElm = f->getNjElm();
//         adrI = f->getAdrI();
//         adrJ = f->getAdrJ();
//         for(int i = 0; i < niElm; ++i){
//           I = adrI[i];
//           if(I < _nInc){
//             for(int j = 0; j < njElm; ++j){
//               J = adrJ[j];
//               if(J < _nInc){
//                 // ierr = MatSetValue(_A, I, J, Ae[njElm*i+j], ADD_VALUES); // TODO : assigner par blocs
//                 ierr = MatSetValue(_A, I, J, Ae[i][j], ADD_VALUES); // TODO : assigner par blocs
//                 CHKERRABORT(PETSC_COMM_WORLD, ierr);
//               }
//             }
//           }
//         }
//       }
//     }
//     ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     // double normMat = 0.0;
//     // ierr = MatNorm(_A, NORM_FROBENIUS, &normMat); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//     // printf("Norme de la matrice : %10.10e\n", normMat);
//   } // if(recomputeMatrix)
//   // viewMatrix();
// }

void feLinearSystemPETSc::assembleMatrices(feSolution *sol){
  PetscErrorCode ierr = 0;
  if(recomputeMatrix){
    // std::cout<<"Computing matrix"<<std::endl;
    std::vector<PetscScalar> values;
    int sizeI, sizeJ;
    std::vector<int> niElm;
    std::vector<int> njElm;
    std::vector<int> adrI;
    std::vector<int> adrJ;

    for(feBilinearForm *f : _formMatrices){
      int nElm = _mesh->getNbElm(f->getCncGeoTag());
      for(int iElm = 0; iElm < nElm; ++iElm){
        f->computeMatrix(_metaNumber, _mesh, sol, iElm); // Matrice elementaire
        double** Ae = f->getAe();
        // Determine assignment indices
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        sizeI = adrI.size(); niElm.reserve(sizeI);
        for(int i = 0; i < sizeI; ++i){ if(adrI[i] < _nInc) niElm.push_back(i); }
        sizeJ = adrJ.size(); njElm.reserve(sizeJ);
        for(int i = 0; i < sizeJ; ++i){ if(adrJ[i] < _nInc) njElm.push_back(i); }
        adrI.erase(std::remove_if( adrI.begin(), adrI.end(), [this](const int& x) { return x >= this->_nInc;  }), adrI.end());
        adrJ.erase(std::remove_if( adrJ.begin(), adrJ.end(), [this](const int& x) { return x >= this->_nInc;  }), adrJ.end());
        // Flatten Ae at relevant indices
        sizeI = adrI.size();
        sizeJ = adrJ.size();
        values.resize(sizeI*sizeJ);
        for(int i = 0; i < sizeI; ++i){
          for(int j = 0; j < sizeJ; ++j){
            values[sizeI*i+j] = Ae[niElm[i]][njElm[j]];
          }
        }
        ierr = MatSetValues(_A, adrI.size(), adrI.data(), adrJ.size(), adrJ.data(), values.data(), ADD_VALUES);
        niElm.clear();
        njElm.clear();
      }
    }
    ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // double normMat = 0.0;
    // ierr = MatNorm(_A, NORM_FROBENIUS, &normMat); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    // printf("Norme de la matrice : %10.10e\n", normMat);
  } // if(recomputeMatrix)
  // MatAXPY(_testA,-1,_A,SAME_NONZERO_PATTERN);
  // viewTestMatrix();
  // MatInfo info;
  // double  mal, nz_a, nz_u, mem, nz_un;

  // MatGetInfo(_A,MAT_LOCAL,&info);
  // mal  = info.mallocs;
  // nz_a = info.nz_allocated;
  // nz_u = info.nz_used;
  // nz_un = info.nz_unneeded;
  // mem = info.memory;

  // std::cout<<"mal = "<<mal<<std::endl;
  // std::cout<<"mem = "<<mal<<std::endl;
  // std::cout<<"nz_a = "<<nz_a<<std::endl;
  // std::cout<<"nz_u = "<<nz_u<<std::endl;
  // std::cout<<"nz_un = "<<nz_un<<std::endl;
  // viewMatrix();
}

// void feLinearSystemPETSc::assembleResiduals(feSolution *sol){
//   PetscErrorCode ierr;
//   PetscInt I;
//   int niElm;
//   std::vector<int> adrI;
//   double normResidual = 0.0;
//   ierr = VecNorm(_res, NORM_2, &normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//   for(feBilinearForm *f : _formResiduals){
//     int nElm = _mesh->getNbElm(f->getCncGeoTag());
//     for(int iElm = 0; iElm < nElm; ++iElm){
//       f->computeResidual(_metaNumber, _mesh, sol, iElm); // Residu elementaire
//       double* Be = f->getBe();
//       niElm = f->getNiElm();
//       adrI = f->getAdrI();
//       for(int i = 0; i < niElm; ++i){
//         I = adrI[i];
//         if(I < _nInc){
//           ierr = VecSetValue(_res, I, Be[i], ADD_VALUES); // TODO : assigner par blocs
//           CHKERRABORT(PETSC_COMM_WORLD, ierr);
//         }
//       }
//     }
//   }
// }

void feLinearSystemPETSc::assembleResiduals(feSolution *sol){
  PetscErrorCode ierr;
  std::vector<PetscScalar> values;
  int sizeI;
  std::vector<int> niElm;
  std::vector<int> adrI;
  double normResidual = 0.0;
  ierr = VecNorm(_res, NORM_2, &normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  for(feBilinearForm *f : _formResiduals){
    int nElm = _mesh->getNbElm(f->getCncGeoTag());
    for(int iElm = 0; iElm < nElm; ++iElm){
      f->computeResidual(_metaNumber, _mesh, sol, iElm); // Residu elementaire
      double* Be = f->getBe();
      // Determine assignment indices
      adrI = f->getAdrI();
      sizeI = adrI.size(); niElm.reserve(sizeI);
      for(int i = 0; i < sizeI; ++i){ if(adrI[i] < _nInc) niElm.push_back(i); }
      adrI.erase(std::remove_if( adrI.begin(), adrI.end(), [this](const int& x) { return x >= this->_nInc;  }), adrI.end());
      // Copy Be into vector
      sizeI = adrI.size();
      values.resize(sizeI);
      for(int i = 0; i < sizeI; ++i){
        values[i] = Be[niElm[i]];
      }
      ierr = VecSetValues(_res, adrI.size(), adrI.data(), values.data(), ADD_VALUES);
      niElm.clear();
    }
  }
  // VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
  // normResidual = 0.0;
  // ierr = VecNorm(_res, NORM_2, &normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // printf("Norme du résidu : %10.10e\n", normResidual);
}

void feLinearSystemPETSc::assemble(feSolution *sol){
  this->assembleMatrices(sol);
  this->assembleResiduals(sol);
}

void feLinearSystemPETSc::solve(double *normDx, double *normResidual, double *normAxb, int *nIter){
  // PetscInt its;
  PetscErrorCode ierr = KSPSolve(ksp,_res,_dx); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  KSPGetIterationNumber(ksp,nIter);
  // PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
  VecSet(_foo, 0.0);
  MatMult(_A,_dx,_foo);
  VecAXPY(_foo, -1.0, _res);
  ierr = VecNorm(_res, NORM_MAX, normResidual); CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm( _dx, NORM_MAX, normDx);       CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = VecNorm( _foo, NORM_MAX, normAxb);       CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // std::cout<<"Norme du residu :        "<<*normResidual<<std::endl;
  // std::cout<<"Norme de la correction : "<<*normDx<<std::endl;
  // std::cout<<"Norme du résidu matriciel Ax-b : "<<normAxb<<std::endl;
  // MatInfo info;
  // double  mal, nz_a, nz_u, mem, nz_un;

  // MatGetInfo(_A,MAT_LOCAL,&info);
  // mal  = info.mallocs;
  // nz_a = info.nz_allocated;
  // nz_u = info.nz_used;
  // nz_un = info.nz_unneeded;
  // mem = info.memory;

  // std::cout<<"mal = "<<mal<<std::endl;
  // std::cout<<"mem = "<<mal<<std::endl;
  // std::cout<<"nz_a = "<<nz_a<<std::endl;
  // std::cout<<"nz_u = "<<nz_u<<std::endl;
  // std::cout<<"nz_un = "<<nz_un<<std::endl;
}

void feLinearSystemPETSc::correctSolution(feSolution *sol){
  std::vector<double> &_sol = sol->getSolutionReference();
  PetscScalar *array;
  VecGetArray(_dx, &array);
  for(int i = 0; i < _nInc; ++i)
    _sol[i] += array[i];
  VecRestoreArray(_dx, &array);
}

void feLinearSystemPETSc::assignResidualToDCResidual(feSolutionContainer *solContainer){
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i)
    solContainer->_fResidual[0][i] = array[i];
  VecRestoreArray(_res, &array);
}

void feLinearSystemPETSc::applyCorrectionToResidual(double coeff, std::vector<double> &d){
  PetscScalar *array;
  VecGetArray(_res, &array);
  for(int i = 0; i < _nInc; ++i)
    array[i] += coeff * d[i];
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