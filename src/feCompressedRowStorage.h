#ifndef _FECOMPRESSEDROWSTORAGE_
#define _FECOMPRESSEDROWSTORAGE_

#include "feNG.h"
#include "feMesh.h"
#include "feNumber.h"
#include "feBilinearForm.h"

class feCompressedRowStorage {
public:

  feCompressedRowStorage(feMetaNumber *metaNumber, feMesh *mesh,
                         std::vector<feBilinearForm *> &formMatrices);
  ~feCompressedRowStorage() {
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
  feInt *liste = NULL;
  bool *ddl_rngcof = NULL;

private:
};

class feCompressedRowStorageMklPardiso : public feCompressedRowStorage {
public:
  feCompressedRowStorageMklPardiso(feMetaNumber *metaNumber, feMesh *mesh,
                                   std::vector<feBilinearForm *> &formMatrices);
  ~feCompressedRowStorageMklPardiso() {
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
  double *allocateMatrix() { return new double[nz]; };
  void freeMatrix(double *Matrix) {
    if(Matrix != NULL) delete[] Matrix;
  };
  void zeroMatrix(double *Matrix) {
    for(feInt i = 0; i < nz; i++) Matrix[i] = 0;
  };
  // ================================================================
  // Somme des coefficients dans la matrice
  // ================================================================
  void matrixAddValues(double *Matrix, feInt nRow, feInt *Row, feInt nColumn, feInt *Column,
                       double **Aij);

protected:
  double *rangee = NULL;
  feInt *irangee = NULL;
  feInt *Ap = NULL;
  feInt *Aj = NULL;

private:
};
#endif