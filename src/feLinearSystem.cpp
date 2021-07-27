#include "feLinearSystem.h"

// void feLinearSystem::viewMatrix(){
//   MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
//   PetscViewer viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,600,600,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView(_A,viewer);
// }

// PetscErrorCode feLinearSystem::initialize(int argc, char** args){
//   PetscErrorCode ierr;
//   PetscMPIInt    size;

//   ierr = PetscInitialize(&argc,&args,(char*) 0,nullptr);
//   if(ierr)
//     return ierr;
//   ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRMPI(ierr);
//   if(size != 1)
//     SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
//   ierr = PetscOptionsGetInt(NULL,NULL,"-n",&_nInc,NULL);    CHKERRQ(ierr);
//   ierr = VecCreate(PETSC_COMM_WORLD,&_dx);                  CHKERRQ(ierr);
//   ierr = PetscObjectSetName((PetscObject) _dx, "Solution"); CHKERRQ(ierr);
//   ierr = VecSetSizes(_dx,PETSC_DECIDE, _nInc);              CHKERRQ(ierr);
//   ierr = VecSetFromOptions(_dx);                            CHKERRQ(ierr);
//   ierr = VecDuplicate(_dx,&_res);                           CHKERRQ(ierr);

//   ierr = MatCreate(PETSC_COMM_WORLD, &_A);              CHKERRQ(ierr);
//   ierr = MatSetSizes(_A,PETSC_DECIDE,PETSC_DECIDE,_nInc,_nInc); CHKERRQ(ierr);
//   ierr = MatSetFromOptions(_A);                         CHKERRQ(ierr);
//   ierr = MatSetUp(_A);                                  CHKERRQ(ierr);

//   ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
//   ierr = KSPSetOperators(ksp,_A,_A);       CHKERRQ(ierr);

//   ierr = KSPGetPC(ksp,&preconditioner);                                         CHKERRQ(ierr);
//   ierr = PCSetType(preconditioner,PCJACOBI);                                       CHKERRQ(ierr);
//   ierr = KSPSetTolerances(ksp,1.e-16,1e-16,PETSC_DEFAULT,100); CHKERRQ(ierr);
//   ierr = KSPSetFromOptions(ksp);                                                CHKERRQ(ierr);
//   return ierr;
// }

// PetscErrorCode feLinearSystem::setToZero(){
//   PetscErrorCode ierr = 0;
//   if(recomputeMatrix){
//     ierr = this->setMatrixToZero();
//     CHKERRQ(ierr);
//   }
//   ierr = this->setResidualToZero();
//   return ierr;
// }

// PetscErrorCode feLinearSystem::setMatrixToZero(){
//   PetscErrorCode ierr = 0;
//   ierr = MatZeroEntries(_A);                      CHKERRQ(ierr);
//   ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//   ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
//   // viewMatrix();
//   return ierr;
// }

// PetscErrorCode feLinearSystem::setResidualToZero(){
//   // viewMatrix();
//   return VecZeroEntries(_res);
// }

// PetscErrorCode feLinearSystem::assembleMatrices(feSolution *sol){
//   PetscErrorCode ierr = 0;
//   if(recomputeMatrix){
//     PetscInt I,J;
//     int niElm, njElm;
//     std::vector<int> adrI, adrJ;

//     for(feBilinearForm *f : _formMatrices){
//       int nElm = _mesh->getNbElm(f->getCncGeoTag());
//       for(int iElm = 0; iElm < nElm; ++iElm){
//         f->computeMatrix(_metaNumber, _mesh, sol, iElm); // Matrice elementaire
//         std::vector<double> Ae = f->getAe();
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
//                 ierr = MatSetValue(_A, I, J, Ae[njElm*i+j], ADD_VALUES); // TODO : assigner par
//                 blocs CHKERRQ(ierr);
//               }
//             }
//           }
//         }
//       }
//     }
//     ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
//   } // if(recomputeMatrix)
//   // viewMatrix();

//   // MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
//   // VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
//   // PetscViewer    viewer;
//   // PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,600,600,&viewer);
//   // PetscObjectSetName((PetscObject)viewer,"Line graph Plot");
//   // PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);

//   // MatOrderingType rtype = MATORDERINGRCM;
//   // IS isrow,iscol;  //row and column permutations
//   // MatGetOrdering(_A,rtype,&isrow,&iscol);
//   // ISView(isrow,PETSC_VIEWER_STDOUT_WORLD);
//   // Mat Aperm;
//   // MatPermute(_A,isrow,iscol,&Aperm);
//   // // VecPermute(_b,colperm,PETSC_FALSE);
//   // MatDestroy(&_A);
//   // _A    = Aperm;
//   // MatView(_A,viewer);
//   // PetscDrawPause();
//   /*
//      Assemble matrix
//   */
//   // PetscInt       i, n = 10, col[3],its;
//   // PetscScalar    value[3];
//   // value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
//   // for (i=1; i<n-1; i++) {
//   //   col[0] = i-1; col[1] = i; col[2] = i+1;
//   //   ierr   = MatSetValues(_A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
//   // }
//   // i    = n - 1; col[0] = n - 2; col[1] = n - 1;
//   // ierr = MatSetValues(_A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
//   // i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
//   // ierr = MatSetValues(_A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
//   // ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   // ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   return ierr;
// }

// PetscErrorCode feLinearSystem::assembleResiduals(feSolution *sol){
//   PetscErrorCode ierr;
//   PetscInt I;
//   int niElm;
//   std::vector<int> adrI;
//   for(feBilinearForm *f : _formResiduals){
//     int nElm = _mesh->getNbElm(f->getCncGeoTag());
//     for(int iElm = 0; iElm < nElm; ++iElm){
//       f->computeResidual(_metaNumber, _mesh, sol, iElm); // Residu elementaire
//       std::vector<double> Be = f->getBe();
//       niElm = f->getNiElm();
//       adrI = f->getAdrI();
//       for(int i = 0; i < niElm; ++i){
//         I = adrI[i];
//         if(I < _nInc){
//           ierr = VecSetValue(_res, I, Be[i], ADD_VALUES); // TODO : assigner par blocs
//           CHKERRQ(ierr);
//         }
//       }
//     }
//   }
//   // VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
//   return ierr;
// }

// PetscErrorCode feLinearSystem::assemble(feSolution *sol){
//   PetscErrorCode ierr;
//   ierr = this->assembleMatrices(sol);
//   ierr = this->assembleResiduals(sol);
//   return ierr;
// }

// PetscErrorCode feLinearSystem::solve(double *normDx, double *normResidual){
//   PetscInt its;
//   PetscErrorCode ierr = KSPSolve(ksp,_res,_dx); CHKERRQ(ierr);
//   // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
//   // KSPGetIterationNumber(ksp,&its);
//   // PetscPrintf(PETSC_COMM_WORLD,"Iterations %D\n",its);
//   ierr = VecNorm(_res, NORM_2, normResidual); CHKERRQ(ierr);
//   ierr = VecNorm( _dx, NORM_2, normDx);       CHKERRQ(ierr);
//   return ierr;
// }

// void feLinearSystem::correctSolution(feSolution *sol){
//   PetscScalar *array;
//   VecGetArray(_dx, &array);
//   for(int i = 0; i < _nInc; ++i)
//     sol->incrementSolAtDOF(i, array[i]);
//   VecRestoreArray(_dx, &array);
// }

// void feLinearSystem::assignResidualToDCResidual(feSolutionContainer *solContainer){
//   PetscScalar *array;
//   VecGetArray(_res, &array);
//   for(int i = 0; i < _nInc; ++i)
//     solContainer->_fResidual[0][i] = array[i];
//   VecRestoreArray(_res, &array);
// }

// void feLinearSystem::printResidual(){
//   // VecView(_res, PETSC_VIEWER_STDOUT_WORLD);
//   PetscScalar *array;
//   VecGetArray(_res, &array);
//   for(int i = 0; i < _nInc; ++i)
//     printf("%12.16f\n", array[i]);
//   VecRestoreArray(_res, &array);
// }

// PetscErrorCode feLinearSystem::check(){
//   VecView(_dx,PETSC_VIEWER_STDOUT_WORLD);

//   PetscReal normRes, normDx;
//   PetscErrorCode ierr = VecNorm(_res,NORM_2,&normRes); CHKERRQ(ierr);
//   ierr = VecNorm( _dx,NORM_2,&normDx); CHKERRQ(ierr);
//   std::cout<<"Norme res = "<<normRes<<std::endl;
//   std::cout<<"Norme dx = "<<normDx<<std::endl;
//   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//   // ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
//   // ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
//   // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
//   // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations
//   %D\n",(double)norm,its);CHKERRQ(ierr); return ierr;
// }

// PetscErrorCode feLinearSystem::finalize(){
//   PetscErrorCode ierr;
//   ierr = KSPDestroy(&ksp);  CHKERRQ(ierr);
//   ierr = MatDestroy(&_A);   CHKERRQ(ierr);
//   ierr = VecDestroy(&_res); CHKERRQ(ierr);
//   ierr = VecDestroy(&_dx);  CHKERRQ(ierr);
//   // ierr = PetscFinalize();
//   return ierr;
// }