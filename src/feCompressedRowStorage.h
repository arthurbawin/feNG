#ifndef _FECOMPRESSEDROWSTORAGE_
#define _FECOMPRESSEDROWSTORAGE_

#include "feNG.h"
#include "feBilinearForm.h"

//
// A simplified class to determine the sparsity pattern (mask)
// of the FE matrix.
// Given to the linear solver to allocate the linear system.
//
class feEZCompressedRowStorage
{
protected:
  // The total number of nonzero entries
  feInt _num_nnz;
  // The number of nonzero for each row (for PETSc)
  std::vector<feInt> nnz;
  // The nonzero entries for each row
  // Using a vector and then sort/unique seems faster than a set
  std::vector<std::vector<feInt> > nnzPerRow;
  // std::vector<std::set<feInt> > nnzPerRow;

  // For parallel PETSc preallocation:
  // The number of nonzero on the diagonal block for each MPI process,
  // and a flattened version to scatter on the procs
  std::vector<std::vector<feInt>> d_nnz;
  std::vector<feInt> d_nnz_flattened;
  // The number of nonzero on the off-diagonal blocks (counted as a whole) for each MPI process
  std::vector<std::vector<feInt>> o_nnz;
  std::vector<feInt> o_nnz_flattened;

  // From Pardiso's doc:
  // ia_Pardiso[i] (i < n) points to the first column index of row i in the array ja_Pardiso
  std::vector<feInt> ia_Pardiso;

  // The (continuous) vector of column indices for each row
  // Columns are given in increasing order and numbered starting at 0
  std::vector<feInt> ja_Pardiso;

public:
  feEZCompressedRowStorage(int numUnknowns, std::vector<feBilinearForm *> &formMatrices,
                           int numMatrixForms, int *ownedUpperBounds = nullptr);
  // ~feEZCompressedRowStorage() = default;

  feInt getNumNNZ() { return _num_nnz; };
  feInt *getNnz() { return nnz.data(); };
  feInt getNnzAt(int i) { return nnz[i]; };

  std::vector<feInt> getDiagNNZForAllMPIProcs() { return d_nnz_flattened; };
  std::vector<feInt> getOffDiagNNZForAllMPIProcs() { return o_nnz_flattened; };

#if defined(HAVE_MKL)
  void get_ia_Pardiso(PardisoInt **array) {
    for(size_t i = 0; i < ia_Pardiso.size(); i++)
      (*array)[i] = (PardisoInt) ia_Pardiso[i];
  }
  void get_ja_Pardiso(PardisoInt **array) { 
  for(size_t i = 0; i < ja_Pardiso.size(); i++)
      (*array)[i] = (PardisoInt) ja_Pardiso[i];
  }
#endif
};

#endif