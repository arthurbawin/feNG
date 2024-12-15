#include "feCompressedRowStorage.h"
#include <mutex>

#if defined(HAVE_PETSC)
#include "petscksp.h"
#endif

// Determine the sparsity pattern of the global FE matrix.
// This is done by doing a "dry" assembly of all the bilinear forms
// associated to a local (elementwise) matrix.
// The sparsity pattern is required by Pardiso,
// and considerably speeds up the PETSc resolution by preallocating the matrix.
// For a parallel (MPI) PETSc matrix, the number of nonzero
// on the diagonal and off-diagonal blocks are also required.
feEZCompressedRowStorage::feEZCompressedRowStorage(int numUnknowns,
                                                   std::vector<feBilinearForm *> &formMatrices,
                                                   int numMatrixForms,
                                                   int *ownedUpperBounds)
{
  nnz.resize(numUnknowns, 0.);
  nnzPerRow.resize(numUnknowns);
  ia_Pardiso.resize(numUnknowns + 1);

  for(int i = 0; i < numUnknowns; ++i) {
    nnzPerRow[i].clear();
    // To add a zero on diagonal elements
    // nnzPerRow[i].insert(i);

    // Difficult to reserve the right vector capacity,
    // the number of indices pushed_back will be proportional
    // to the number of forms, 
    nnzPerRow[i].push_back(i);
  }

  for(int iForm = 0; iForm < numMatrixForms; ++iForm)
  {
    // feBilinearForm *f = formMatrices[iForm];
    feCncGeo *cnc = formMatrices[iForm]->getCncGeo();
    int numColors = cnc->getNbColor();
    const std::vector<int> &numElemPerColor = cnc->getNbElmPerColor();
    const std::vector<std::vector<int> > &listElmPerColor = cnc->getListElmPerColor();
    int numElementsInColor;
    std::vector<int> listElmC;

    for(int iColor = 0; iColor < numColors; ++iColor)
    {
      numElementsInColor = numElemPerColor[iColor];
      listElmC = listElmPerColor[iColor];

      #if defined(HAVE_OMP)
      #pragma omp parallel
      #endif
      {
        int elm, nI, nJ;

        #if defined(HAVE_OMP)
        #pragma omp for
        #endif
        for(int iElm = 0; iElm < numElementsInColor; ++iElm)
        {
          #if defined(HAVE_OMP)
          int numForm = iForm + omp_get_thread_num() * numMatrixForms;
          feBilinearForm *f = formMatrices[numForm];
          #else
          feBilinearForm *f = formMatrices[iForm];
          #endif

          elm = listElmC[iElm];

          f->initializeAddressingVectors(elm);

          // Determine global assignment indices
          std::vector<feInt> &adrI = f->getAdrI();
          std::vector<feInt> &adrJ = f->getAdrJ();
          nI = adrI.size();
          nJ = adrJ.size();

          for(int i = 0; i < nI; ++i) {
            for(int j = 0; j < nJ; ++j) {
              if(adrI[i] < numUnknowns && adrJ[j] < numUnknowns)
              {
                // If using a vector of sets
                // nnzPerRow[adrI[i]].insert(adrJ[j]);

                // This should be thread-safe as we loop
                // over elements of the same color
                nnzPerRow[adrI[i]].push_back(adrJ[j]);
              }
            }
          }
        }
      }
    }
  }

  // Sort and unique
  #if defined(HAVE_OMP)
  #pragma omp parallel for
  #endif
  for(int i = 0; i < numUnknowns; ++i) {
    std::sort( nnzPerRow[i].begin(), nnzPerRow[i].end() );
    nnzPerRow[i].erase( std::unique( nnzPerRow[i].begin(), nnzPerRow[i].end() ), nnzPerRow[i].end() );
  }

  _num_nnz = 0;
  for(int i = 0; i < numUnknowns; ++i) {
    nnz[i] = fmax(1, nnzPerRow[i].size());
    _num_nnz += nnz[i];
  }

  // Copy nnzPerRow size_to continuous vector
  ja_Pardiso.resize(_num_nnz);
  feInt cnt = 0;
  for(int i = 0; i < numUnknowns; ++i) {
    ia_Pardiso[i] = cnt;
    for(auto &val : nnzPerRow[i]) {
      ja_Pardiso[cnt++] = val;
    }
  }
  ia_Pardiso[numUnknowns] = _num_nnz;

  // Also determine diagonal and off-diagonal nnz for parallel PETSc matrix
  if(ownedUpperBounds != nullptr) {
#if defined(HAVE_PETSC)
    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    d_nnz.resize(size);
    o_nnz.resize(size);

    // for(int i = 0; i < size; ++i) {
    //   feInfo("Bound = %d", ownedUpperBounds[i]);
    // }

    int iBound = 0;
    int currentProc = 0, lowerBound = 0, upperBound = ownedUpperBounds[iBound];
    d_nnz_flattened.clear();
    o_nnz_flattened.clear();

    // Loop over the rows
    for(int i = 0; i < numUnknowns; ++i) {

      // Go to the structure on the next proc if necessary
      if(i >= upperBound) {
        iBound++;
        lowerBound = upperBound;
        upperBound = ownedUpperBounds[iBound];
        currentProc++;
      }

      // feInfo("Row %d - Owned by proc %d - num total nnz = %d ", i, currentProc, nnzPerRow[i].size());

      int nnzOnRowOfDiagBlock = 0, nnzOnRowOfOffDiagBlock = 0;
      for(auto val : nnzPerRow[i]) {

        // feInfo("%d", val);

        if(lowerBound <= val && val < upperBound) {
          // Nonzero entry is in the diagonal proc owned by this proc
          nnzOnRowOfDiagBlock++;
        } else {
          nnzOnRowOfOffDiagBlock++;
        }

      }

      // feInfo("%d nnz on diag block - %d on off-diag", nnzOnRowOfDiagBlock, nnzOnRowOfOffDiagBlock);
      d_nnz[currentProc].push_back(nnzOnRowOfDiagBlock);
      o_nnz[currentProc].push_back(nnzOnRowOfOffDiagBlock);
      d_nnz_flattened.push_back(nnzOnRowOfDiagBlock);
      o_nnz_flattened.push_back(nnzOnRowOfOffDiagBlock);
    }

#else
    d_nnz.resize(1);
    o_nnz.resize(1);
#endif
  }
}