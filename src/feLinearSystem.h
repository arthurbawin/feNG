#ifndef _FELINEARSYSTEM_
#define _FELINEARSYSTEM_

#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"

#include "petscksp.h"

class feLinearSystem{

protected:
  std::vector<feBilinearForm*> &_formMatrices;
  std::vector<feBilinearForm*> &_formResiduals;
  feMetaNumber *_metaNumber;
  feMesh *_mesh;

  PetscInt _nInc;
  PetscInt _nDofs;
  Mat _A;
  Vec _res;
  Vec _dx;

  KSP ksp;
  PC preconditioner;

  bool recomputeMatrix;
public:
	feLinearSystem(std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
    feMetaNumber *metaNumber, feMesh *mesh)  : _formMatrices(formMatrices), _formResiduals(formResiduals), 
  _metaNumber(metaNumber), _mesh(mesh), _nInc(metaNumber->getNbUnknowns()), _nDofs(metaNumber->getNbDOFs()),
  recomputeMatrix(false)
  {
    
  };
	~feLinearSystem() {}

  bool getRecomputeStatus(){ return recomputeMatrix; }
  void setRecomputeStatus(bool status){ recomputeMatrix = status; }

  void viewMatrix();

  PetscErrorCode initialize(int argc, char** args);
  PetscErrorCode finalize();

  PetscErrorCode setToZero();

  PetscErrorCode assembleMatrices(feSolution *sol);
  PetscErrorCode assembleResiduals(feSolution *sol);
  PetscErrorCode assemble(feSolution *sol);

  PetscErrorCode solve(double *normDx, double *normResidual);
  void correctSolution(feSolution *sol);
  PetscErrorCode check();
};

#endif