#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"

#include "petscksp.h"

class feLinearSystem{

protected:
  std::vector<feBilinearForm*> &_forms;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  PetscInt _nInc;
  PetscInt _nDofs;
  Mat _A;
  Vec _res;
  Vec _dx;

  KSP ksp;
  PC preconditioner;

public:
	feLinearSystem(std::vector<feBilinearForm*> &forms, feMetaNumber *metaNumber, feMesh *mesh)
    : _forms(forms), _metaNumber(metaNumber), _mesh(mesh), _nInc(metaNumber->getNbUnknowns()), 
    _nDofs(metaNumber->getNbDOFs())
  {
    
  };
	~feLinearSystem() {}

  PetscErrorCode initialize(int argc, char** args){
    PetscErrorCode ierr;
    PetscMPIInt    size;

    ierr = PetscInitialize(&argc,&args,(char*) 0,nullptr);
    if(ierr)
      return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRMPI(ierr);
    if(size != 1)
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_WRONG_MPI_SIZE,"This is a uniprocessor example only!");
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&_nInc,NULL);    CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&_dx);                  CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) _dx, "Solution"); CHKERRQ(ierr);
    ierr = VecSetSizes(_dx,PETSC_DECIDE, _nInc);              CHKERRQ(ierr);
    ierr = VecSetFromOptions(_dx);                            CHKERRQ(ierr);
    ierr = VecDuplicate(_dx,&_res);                           CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &_A);              CHKERRQ(ierr);
    ierr = MatSetSizes(_A,PETSC_DECIDE,PETSC_DECIDE,_nInc,_nInc); CHKERRQ(ierr);
    ierr = MatSetFromOptions(_A);                         CHKERRQ(ierr);
    ierr = MatSetUp(_A);                                  CHKERRQ(ierr);

    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,_A,_A);       CHKERRQ(ierr);

    ierr = KSPGetPC(ksp,&preconditioner);                                         CHKERRQ(ierr);
    ierr = PCSetType(preconditioner,PCJACOBI);                                    CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp);                                                CHKERRQ(ierr);
  }

  PetscErrorCode assemble(feSolution *sol){
    PetscErrorCode ierr;
    PetscInt I,J;
    int niElm, njElm;
    std::vector<int> adrI, adrJ;

    for(feBilinearForm *f : _forms){
      int nElm = _mesh->getNbElm(f->getCncGeoTag());
      for(int iElm = 0; iElm < nElm; ++iElm){
        // Matrice et residu elementaires
        f->computeMatrix(_metaNumber, _mesh, sol, iElm);
        f->computeRHS(_metaNumber, _mesh, sol, iElm);
        // Assemblage
        std::vector<double> Ae = f->getAe(), Be = f->getBe();
        niElm = f->getNiElm();
        njElm = f->getNjElm();
        adrI = f->getAdrI();
        adrJ = f->getAdrJ();
        for(int i = 0; i < niElm; ++i){
          I = adrI[i];
          if(I < _nInc){
            ierr = VecSetValue(_res, I, Be[i], ADD_VALUES); // TODO : assigner par blocs
            CHKERRQ(ierr);
            for(int j = 0; j < njElm; ++j){
              J = adrJ[j];
              if(J < _nInc){
                ierr = MatSetValue(_A, I, J, Ae[njElm*i+j], ADD_VALUES); // TODO : assigner par blocs
                CHKERRQ(ierr);
              }
            }
          }
        }
      }
    }
    ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);   CHKERRQ(ierr);
    MatView(_A,PETSC_VIEWER_STDOUT_WORLD);
    VecView(_res,PETSC_VIEWER_STDOUT_WORLD);
    /*
       Assemble matrix
    */
    // PetscInt       i, n = 10, col[3],its;
    // PetscScalar    value[3];
    // value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
    // for (i=1; i<n-1; i++) {
    //   col[0] = i-1; col[1] = i; col[2] = i+1;
    //   ierr   = MatSetValues(_A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    // }
    // i    = n - 1; col[0] = n - 2; col[1] = n - 1;
    // ierr = MatSetValues(_A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    // i    = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
    // ierr = MatSetValues(_A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    // ierr = MatAssemblyBegin(_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // ierr = MatAssemblyEnd(_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  PetscErrorCode solve(){
    PetscErrorCode ierr = KSPSolve(ksp,_res,_dx); CHKERRQ(ierr);
  }

  PetscErrorCode check(){
    VecView(_dx,PETSC_VIEWER_STDOUT_WORLD);

    PetscReal normRes, normDx;
    PetscErrorCode ierr = VecNorm(_res,NORM_2,&normRes); CHKERRQ(ierr);
    ierr = VecNorm( _dx,NORM_2,&normDx); CHKERRQ(ierr);
    std::cout<<"Norme res = "<<normRes<<std::endl;
    std::cout<<"Norme dx = "<<normDx<<std::endl;
    // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    // ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
  }
  
  PetscErrorCode finalize(){
    PetscErrorCode ierr;
    ierr = VecDestroy(&_dx);  CHKERRQ(ierr);
    ierr = VecDestroy(&_res); CHKERRQ(ierr); 
    ierr = MatDestroy(&_A);   CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);  CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
  }
};

#endif