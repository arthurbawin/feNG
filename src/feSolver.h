#ifndef _FESOLVER_
#define _FESOLVER_

#include "feMesh.h"
#include "feLinearSystem.h"
#include "feBilinearForm.h"
#include "feNumber.h"
#include "feNorm.h"

typedef struct feTol{
  double tolDx;
  double tolResidual;
  double maxIter;
} feTolerances;

double solveStationary(feTolerances tol, feMetaNumber *metaNumber, feLinearSystem *linearSystem, 
  std::vector<feBilinearForm*> &formMatrices, std::vector<feBilinearForm*> &formResiduals,
  feSolution *sol, feNorm *norm, feMesh *mesh);

// class feSolution{

// protected:

// public:
// 	feSolution()
// 	{

// 	};
// 	~feSolution() {
//   }


// };

#endif