#ifndef _FELINEARSYSTEMMKLPARDISO_
#define _FELINEARSYSTEMMKLPARDISO_

#include "feNumber.h"
#include "feMesh.h"
#include "feBilinearForm.h"
#include "feSolution.h"
#include "feLinearSystem.h"
#include "feCompressedRowStorage.h"
#include "feSolutionContainer.h"

#if defined(HAVE_MKL)

class feLinearSystemMklPardiso : public feLinearSystem
{
public:
  feLinearSystemMklPardiso(std::vector<feBilinearForm *> bilinearForms, feMetaNumber *metaNumber,
                           feMesh *mesh);
  virtual ~feLinearSystemMklPardiso();

  // ====================================================================
  // Pour modifier les paramètres de Pardiso
  // pivot : valeur de l'exposant du nombre réel 10^(-pivot)
  //         valeur par défaut pivot := 13
  // ====================================================================
  void setPivot(int pivot);
  void setPardisoMsglvlHigh() { MSGLVL = 1; };
  void setPardisoMsglvlLow() { MSGLVL = 0; };
  void print_matrix();
  bool getSymbolicFactorizationStatus() { return symbolicFactorization; }
  void setSymbolicFactorizationStatus(bool status)
  {
    symbolicFactorization = status;
    if(symbolicFactorization) recomputeMatrix = true;
  }
  void toggleSymbolicFactorizationStatus()
  {
    symbolicFactorization = !symbolicFactorization;
    if(symbolicFactorization) recomputeMatrix = true;
  }
  // ====================================================================
  // Méthodes pour controler le calcul de la matrice
  // solve              : résoudre le système
  // setTozero          : initialise la matrice et le résidu à zéro
  // setMatrixToZero    : initialise la matrice à zéro
  // setResidualToZero  : initialise le résidu à zéro
  // assemble           : somme de la matrice et du résidu
  // assembleMatrix     : somme de la matrice
  // assembleResidual   : somme du résidu
  // correctSolution    : somme de la correction dans le vecteur solution
  //                      ou dans la structure solution
  // ====================================================================
  virtual void solve(double *normDx, double *normResidual, double *normAxb,
                     int *nIter); // MATLAB : RESOUDRE
  virtual void setToZero(); // MATLAB : INITIALISE
  virtual void setMatrixToZero(); // MATLAB : INITIALISE_LA_MATRICE
  virtual void setResidualToZero(); // MATLAB : INITIALISE_LE_RESIDU
  virtual void assembleMatrices(feSolution *sol); // MATLAB : SOMME_LA_MATRICE
  virtual void assembleResiduals(feSolution *sol); // MATLAB : SOMME_LE RESIDU
  virtual void assemble(feSolution *sol); // MATLAB : SOMME MATRICE ET RESIDU
  virtual void assignResidualToDCResidual(feSolutionContainer *solContainer);
  virtual void applyCorrectionToResidual(double coeff, std::vector<double> &d);
  virtual void correctSolution(feSolution *sol); // MATLAB : ------
  virtual void correctSolution(double *sol); // MATLAB : ------
protected:
  //=====================================================================
  // Structure Creuse CSR de MKL
  //=====================================================================

  feCompressedRowStorageMklPardiso *crsMklPardiso;
  feMKLPardisoInt matrixOrder;
  double *du;
  double *residu;
  feInt iparm12 = 1; // Modification de IPARM[12]
  bool symbolicFactorization = true;

  //=====================================================================
  // STRUCTURES DE DONNEES INTERNES DE PARDISO-MKL
  //=====================================================================
  // =PARDISOINIT===================
  void *PT[64];
  feMKLPardisoInt MYTPE;
  feMKLPardisoInt IPARM[64];
  double DPARM[64];

  // =PARDISO=======================
  feMKLPardisoInt MAXFCT;
  feMKLPardisoInt MNUM;
  feMKLPardisoInt MTYPE;
  feMKLPardisoInt SOLVER;
  feMKLPardisoInt PHASE;
  feMKLPardisoInt N;
  feMKLPardisoInt NRHS;
  feMKLPardisoInt ERROR;
  feMKLPardisoInt MSGLVL;
  feMKLPardisoInt IDUM;
  double DDUM;
  feMKLPardisoInt IPIVOT;

  feMKLPardisoInt nz;
  feMKLPardisoInt *Ap; // dimension ordre+1, IA
  feMKLPardisoInt *Aj; // dimension nz, JA
  double *Ax; // dimension nz

private:
  //=====================================================================
  // Les méthodes de Pardiso
  //=====================================================================
  void mklSymbolicFactorization(void);
  void mklFactorization(void);
  void mklSolve(void);
};

#endif // HAVE_MKL

#endif