#include "feCompressedRowStorage.h"
#include <mutex>

feEZCompressedRowStorage::feEZCompressedRowStorage(int numUnknowns,
                                                   std::vector<feBilinearForm *> &formMatrices,
                                                   int numMatrixForms)
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

  for(size_t iForm = 0; iForm < numMatrixForms; ++iForm)
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
  #pragma omp parallel for
  for(size_t i = 0; i < numUnknowns; ++i) {
    std::sort( nnzPerRow[i].begin(), nnzPerRow[i].end() );
    nnzPerRow[i].erase( std::unique( nnzPerRow[i].begin(), nnzPerRow[i].end() ), nnzPerRow[i].end() );
  }

  _num_nnz = 0;
  for(size_t i = 0; i < numUnknowns; ++i) {
    nnz[i] = fmax(1, nnzPerRow[i].size());
    _num_nnz += nnz[i];
  }

  // Copy nnzPerRow size_to continuous vector
  ja_Pardiso.resize(_num_nnz);
  feInt cnt = 0;
  for(int i = 0; i < numUnknowns; ++i) {
    ia_Pardiso[i] = cnt;
    for(auto &val : nnzPerRow[i])
      ja_Pardiso[cnt++] = val;
  }
  ia_Pardiso[numUnknowns] = _num_nnz;
}


//
// Everything below is deprecated
//

// ====================================================================
// Constructeur de la classe de base CSR
// ====================================================================

feCompressedRowStorage::feCompressedRowStorage(int numUnknowns,
                                               std::vector<feBilinearForm *> &formMatrices,
                                               int numMatrixForms)
{
  ordre = (feInt)numUnknowns;
  nnz = new feInt[ordre];

  ddlNumberOfElements = new feInt[ordre]; // à initialiser à zéro
  for(feInt i = 0; i < ordre; i++) ddlNumberOfElements[i] = 0;

  for(feBilinearForm *f : formMatrices) {
    feInt nbElems = f->getCncGeo()->getNumElements();

    for(feInt iElm = 0; iElm < nbElems; iElm++) {
      f->initializeAddressingVectors(iElm);
      feInt NBRI = f->getLocalMatrixM();
      std::vector<feInt> VADI = f->getAdrI(); // &VADI ????

      for(feInt i = 0; i < NBRI; i++) {
        feInt I = (feInt)VADI[i];
        if(I < ordre) ddlNumberOfElements[I]++;
      }
    }
  }

  // for(int i = 0; i < ordre; ++i){
  //   feInfo("ddl[%d] = %d", i, ddlNumberOfElements[i]);
  // }

  feInt totalNumberOfElements = 0;
  for(feInt i = 0; i < ordre; i++) totalNumberOfElements += ddlNumberOfElements[i];
  // ================================================================
  // IDENTIFIER LES ELEMENTS ET LES FORMES BILINEAIES D'UN
  // DEGRE DE LIBERTE -- PREPARARION DE L'ESPACE DE TRAVAIL
  // ================================================================
  ddlBiLinearForms = new feInt *[ordre];
  ddlElms = new feInt *[ordre];

  feInt *ptr = new feInt[totalNumberOfElements];
  for(feInt i = 0; i < ordre; i++) {
    ddlBiLinearForms[i] = ptr;
    ptr += ddlNumberOfElements[i];
  }

  ptr = new feInt[totalNumberOfElements];
  for(feInt i = 0; i < ordre; i++) {
    ddlElms[i] = ptr;
    ptr += ddlNumberOfElements[i];
  }

  // ================================================================
  // RECUPERER LES NUMEROS DES FORMES BILINEAIRES ET DES ELEMENTS
  // ================================================================
  for(feInt i = 0; i < ordre; i++) ddlNumberOfElements[i] = 0;
  for(int eq = 0; eq < numMatrixForms; eq++) {
    feBilinearForm *equelm = formMatrices[eq];
    feInt nbElems = equelm->getCncGeo()->getNumElements();

    for(feInt el = 0; el < nbElems; el++) {
      equelm->initializeAddressingVectors(el);
      feInt NBRI = equelm->getLocalMatrixM();
      std::vector<feInt> VADI = equelm->getAdrI();

      for(feInt i = 0; i < NBRI; i++) {
        feInt I = VADI[i];
        if(I < ordre) {
          feInt k = ddlNumberOfElements[I];
          ddlBiLinearForms[I][k] = eq;
          ddlElms[I][k] = el;
          ddlNumberOfElements[I]++;
        }
      }
    }
  }
  // ================================================================
  // COMPTER LES COEFFICIENTS D'UNE RANGEE
  // ================================================================
  // feInt* liste        = new feInt [ordre];
  // bool*  ddl_rngcof   = new bool  [ordre];
  // feInt* ddl_nbrcof   = new feInt [ordre];
  liste = new feInt[ordre];
  ddl_rngcof = new bool[ordre];
  feInt *ddl_nbrcof = new feInt[ordre];
  feInt eq = 0;
  feInt el = 0;

  for(feInt i = 0; i < ordre; i++) {
    ddl_rngcof[i] = false;
    liste[i] = 0;
  }
  for(feInt I = 0; I < ordre; I++) {
    feInt nbrelm = ddlNumberOfElements[I];
    feInt nbrcof = 0;

    // ================================================================
    // Pour toujours avoir un coefficient sur la diagonale
    // 18 fevrier 2011 mod -- AG
    // ================================================================
    ddl_rngcof[I] = true;
    liste[nbrcof] = I;
    nbrcof++;

    for(feInt k = 0; k < nbrelm; k++) {
      eq = ddlBiLinearForms[I][k];
      el = ddlElms[I][k];

      feBilinearForm *equelm = formMatrices[eq];
      // feInt cncGeoTag = equelm->getCncGeoTag();
      equelm->initializeAddressingVectors(el);
      feInt NBRJ = equelm->getLocalMatrixM();
      std::vector<feInt> VADJ = equelm->getAdrJ();

      for(feInt j = 0; j < NBRJ; j++) {
        feInt J = VADJ[j];
        if(J < ordre) {
          if(!ddl_rngcof[J]) {
            ddl_rngcof[J] = true;
            liste[nbrcof] = J;
            nbrcof++;
          }
        }
      }
    }

    ddl_nbrcof[I] = nbrcof;
    for(feInt i = 0; i < nbrcof; i++) ddl_rngcof[liste[i]] = false;
  }

  // ================================================================
  // ON CONSERVE NNZ
  // ================================================================

  for(feInt i = 0; i < ordre; i++) nnz[i] = ddl_nbrcof[i];
  nz = 0;
  for(feInt i = 0; i < ordre; i++) nz += nnz[i];

  // ================================================================
  // LIBERER L'ESPACE DE TRAVAIL TEMPORAIRE
  // ================================================================
  // delete [] ddlNumberOfElements;

  // delete [] ddlBiLinearForms[0];
  // delete [] ddlBiLinearForms;
  // delete [] ddlElms[0];
  // delete [] ddlElms;

  // delete [] liste;
  // delete [] ddl_rngcof;
  delete[] ddl_nbrcof;
}

// ====================================================================
// Constructeur de la classe dérivée CSR MKL PARDISO
// ====================================================================

//========================================================================================================
// fonction pour comprarer deux entiers, pour le tri
// d'une liste d'entiers par la fonction qsort
//========================================================================================================
int compint(const void *a, const void *b)
{
  const feInt *aa = (const feInt *)a;
  const feInt *bb = (const feInt *)b;
  int code = 0;
  feInt val = (*aa - *bb);
  if(val < 0)
    code = -1;
  else if(val == 0)
    code = 0;
  else
    code = 1;
  return code;
}

feCompressedRowStorageMklPardiso::feCompressedRowStorageMklPardiso(
  int numUnknowns, std::vector<feBilinearForm *> &formMatrices, int numMatrixForms)
  : feCompressedRowStorage(numUnknowns, formMatrices, numMatrixForms)
{
  // ================================================================
  // ALLOUER LA STRUCTURE CSR DE MKL PARDISO
  // ================================================================
  irangee = new feInt[ordre];
  rangee = new double[ordre];
  Ap = new feInt[ordre + 1];
  Aj = new feInt[nz]; // EFMKLPARDISOint
  // Ax      = new double[nz];
  // ================================================================
  // INITIALISATION
  // ================================================================
  for(feInt i = 0; i < ordre; i++) rangee[i] = 0;
  for(feInt i = 0; i < ordre + 1; i++) Ap[i] = 0;
  for(feInt i = 0; i < ordre; i++) Aj[i] = 0;

  // =============================================================
  // GENERER LA STRUCTURE Ap
  // =============================================================
  if(Ap != NULL) {
    Ap[0] = 0;
    for(feInt i = 0; i < ordre; i++) {
      // feInfo("nnz[%d] = %d", i, nnz[i]);
      Ap[i + 1] = Ap[i] + (feInt)nnz[i]; // EFMKLPARDISOint
    }
  }

  // ================================================================
  // COMPTER LES COEFFICIENTS D'UNE RANGEE
  // ON COMPTE DE NOUVEAU POUR CONSERVER LES NUMÉROS DES COLONNES
  // DANS LA STRUCTURE Aj
  // ================================================================
  feInt eq = 0;
  feInt el = 0;

  for(feInt i = 0; i < ordre; i++) {
    ddl_rngcof[i] = false;
    liste[i] = 0;
  }
  for(feInt I = 0; I < ordre; I++) {
    // Un pinteur sur le début de la rangée
    feInt *ptrAj = Aj + Ap[I];

    feInt nbrelm = ddlNumberOfElements[I];
    feInt nbrcof = 0;

    // ================================================================
    // Pour stoujours avoir un coefficient sur la diagonale
    // 18 fevrier 2011 mod -- AG
    // ================================================================
    ddl_rngcof[I] = true;
    liste[nbrcof] = I;
    nbrcof++;

    for(feInt k = 0; k < nbrelm; k++) {
      eq = ddlBiLinearForms[I][k];
      el = ddlElms[I][k];

      feBilinearForm *equelm = formMatrices[eq];
      // feInt cncGeoTag = equelm->getCncGeoTag();
      equelm->initializeAddressingVectors(el);
      feInt NBRJ = equelm->getLocalMatrixM();
      std::vector<feInt> VADJ = equelm->getAdrJ();

      for(feInt j = 0; j < NBRJ; j++) {
        feInt J = VADJ[j];
        if(J < ordre) {
          if(!ddl_rngcof[J]) {
            ddl_rngcof[J] = true;
            liste[nbrcof] = J;
            nbrcof++;
          }
        }
      }
    }

    for(feInt i = 0; i < nbrcof; i++) ddl_rngcof[liste[i]] = false;

    // Il faut conserver la liste des colonnes
    qsort(liste, nbrcof, sizeof(feInt), compint);

    for(feInt i = 0; i < nbrcof; i++) {
      ddl_rngcof[liste[i]] = false;
      *ptrAj = (feInt)liste[i];
      ptrAj++;
    }
  }

  // ================================================================
  // PARDISO
  // ================================================================
  // ================================================================
  // CONVERSION BASE 0 C A BASE 1 FORTRAN
  // ================================================================
  for(feInt i = 0; i < nz; i++) Aj[i] += 1;
  for(feInt i = 0; i < ordre + 1; i++) Ap[i] += 1;
}

void feCompressedRowStorageMklPardiso::print_info()
{
  // for(feInt i=0;i<nz;i++) printf(" %ld \n", Aj[i]);
  for(feInt i = 0; i < ordre; i++) {
    printf("rangee (%d) ", i + 1);
    for(feInt j = 0; j < Ap[i + 1] - Ap[i]; j++) printf(" %d ", Aj[Ap[i] - 1 + j]);
    printf("\n");
  }
}

// void feCompressedRowStorageMklPardiso::matrixAddValues(double *Matrix, feInt nRow, feInt *Row,
//                                                        feInt nColumn, feInt *Column, double **Aij)
// {
//   for(feInt i = 0; i < nRow; i++) {
//     feInt I = Row[i];
//     if(I < ordre) {
//       feInt debut = Ap[I] - 1;
//       feInt fin = Ap[I + 1] - 1;
//       feInt ncf = fin - debut;

//       for(feInt j = 0; j < ncf; j++) irangee[Aj[debut + j] - 1] = debut + j;
//       for(feInt j = 0; j < nColumn; j++) {
//         feInt J = Column[j];
//         if(J < ordre) {
//           // printf("irangee %ld Aij %g\n", irangee[J], Aij[i][j]);
//           Matrix[irangee[J]] += Aij[i][j];
//         }
//       }
//     }
//   }
// }