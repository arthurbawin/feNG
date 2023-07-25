#ifndef _FECOMPRESSEDROWSTORAGE_
#define _FECOMPRESSEDROWSTORAGE_

#include "feNG.h"
#include "feBilinearForm.h"

// A simplified class to give to PETSc
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

  // From Pardiso's doc:
  // ia_Pardiso[i] (i < n) points to the first column index of row i in the array ja_Pardiso
  std::vector<feInt> ia_Pardiso;

  // The (continuous) vector of column indices for each row
  // Columns are given in increasing order and numbered starting at 0
  std::vector<feInt> ja_Pardiso;

public:
  feEZCompressedRowStorage(int numUnknowns, std::vector<feBilinearForm *> &formMatrices,
                           int numMatrixForms);
  ~feEZCompressedRowStorage(){};

  feInt getNumNNZ() { return _num_nnz; };
  feInt *getNnz() { return nnz.data(); };
  feInt getNnzAt(int i) { return nnz[i]; };

#if defined(HAVE_MKL)
  void get_ia_Pardiso(PardisoInt **array) {
    for(feInt i = 0; i < ia_Pardiso.size(); i++)
      (*array)[i] = (PardisoInt) ia_Pardiso[i];
  }
  void get_ja_Pardiso(PardisoInt **array) { 
  for(feInt i = 0; i < ja_Pardiso.size(); i++)
      (*array)[i] = (PardisoInt) ja_Pardiso[i];
  }
#endif

  double *allocateMatrixArray() { return new double[_num_nnz]; };
  void freeMatrixArray(double *M)
  {
    if(M != NULL) delete[] M;
  };
  void setMatrixArrayToZero(double *M)
  {
    for(feInt i = 0; i < _num_nnz; i++) M[i] = 0.;
  };
};

// Deprecated and should be removed, used only in the eigensolver
class feCompressedRowStorage
{
public:
  feCompressedRowStorage(int numUnknowns, std::vector<feBilinearForm *> &formMatrices,
                         int numMatrixForms);
  ~feCompressedRowStorage()
  {
    if(nnz != NULL) delete[] nnz;
    // ================================================================
    // Structure temporaire nécessaire à la construction CSR
    // Utile pour les classes dérivées
    // ================================================================
    if(ddlNumberOfElements != NULL) delete[] ddlNumberOfElements;
    if(ddlBiLinearForms != NULL) {
      if(ddlBiLinearForms[0] != NULL) delete[] ddlBiLinearForms[0];
      delete[] ddlBiLinearForms;
    }

    if(ddlElms != NULL) {
      if(ddlElms[0] != NULL) delete[] ddlElms[0];
      delete[] ddlElms;
    }

    if(ddlColor != NULL) {
      if(ddlColor[0] != NULL) delete[] ddlColor[0];
      delete[] ddlColor;
    }

    if(liste != NULL) delete[] liste;
    if(ddl_rngcof != NULL) delete[] ddl_rngcof;
  };

  feInt *getNnz() { return nnz; };
  feInt getNz() { return nz; };
  feInt getMatrixOrder() { return ordre; };

protected:
  // ================================================================
  // ordre ; la taille de la matrice carrée ordre x ordre
  // nnz   ; le nombre de coefficients non nuls par rangée
  // nz    ; le nombre total de coefficients non nuls
  // ================================================================
  feInt ordre = 0;
  feInt *nnz = NULL;
  feInt nz = 0;

  // ================================================================
  // Structure temporaire nécessaire à la construction CSR
  // Elles sont conservées car elles sont utiles pour les classes
  // dérivées
  // ================================================================
  feInt *ddlNumberOfElements = NULL;
  feInt **ddlBiLinearForms = NULL;
  feInt **ddlElms = NULL;
  feInt **ddlColor = NULL;
  feInt *liste = NULL;
  bool *ddl_rngcof = NULL;
};

class feCompressedRowStorageMklPardiso : public feCompressedRowStorage
{
public:
  feCompressedRowStorageMklPardiso(int numUnknowns, std::vector<feBilinearForm *> &formMatrices,
                                   int numMatrixForms);
  ~feCompressedRowStorageMklPardiso()
  {
    if(rangee != NULL) delete[] rangee;
    if(irangee != NULL) delete[] irangee;
    if(Ap != NULL) delete[] Ap;
    if(Aj != NULL) delete[] Aj;
  };
  void print_info();

  // ================================================================
  // Interface avec MKL PARDISO
  // ================================================================
  feInt *getAp() { return Ap; };
  feInt *getAj() { return Aj; };
  // double *getRangee() { return rangee; };
  // feInt *getIrangee() { return irangee; };
  double *allocateMatrix() { return new double[nz]; };
  void freeMatrix(double *Matrix)
  {
    if(Matrix != NULL) delete[] Matrix;
  };
  void zeroMatrix(double *Matrix)
  {
    for(feInt i = 0; i < nz; i++) Matrix[i] = 0;
  };
  // ================================================================
  // Somme des coefficients dans la matrice
  // ================================================================
  // void matrixAddValues(double *Matrix, feInt nRow, feInt *Row, feInt nColumn, feInt *Column,
  //                      double **Aij);

protected:
  double *rangee = NULL;
  feInt *irangee = NULL;
  feInt *Ap = NULL;
  feInt *Aj = NULL;

private:
};
#endif