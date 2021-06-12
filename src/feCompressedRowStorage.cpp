#include "feCompressedRowStorage.h"

// ====================================================================
// Constructeur de la classe de base CSR
// ====================================================================
feCompressedRowStorage::feCompressedRowStorage(feMetaNumber *metaNumber, feMesh *mesh, std::vector<feBilinearForm*> &formMatrices){
  ordre = (feInt) metaNumber->getNbUnknowns();
  nnz   = new feInt [ordre];

  ddlNumberOfElements = new feInt [ordre];  // à initialiser à zéro
  for(feInt i=0;i<ordre;i++) ddlNumberOfElements[i]=0;

  // ================================================================
  // COMPTER LES ELEMENTS D'UN DEGRE DE LIBERTE
  // POUR CONSTRUIRE L'ESPACE TEMPORAIRE 
  // ================================================================
  feInt NumberOfBilinearForms = formMatrices.size();

  for(feInt eq=0;eq<NumberOfBilinearForms;eq++){
    feBilinearForm* equelm = formMatrices[eq];
    feInt cncGeoTag = equelm->getCncGeoTag();
    feInt nbElems   = mesh->getNbElm(cncGeoTag);
    for(feInt e=0;e<nbElems;e++) {
      equelm->initialize_vadij_only(metaNumber, e);
      feInt             NBRI = equelm->getNiElm();
      std::vector<int>  VADI = equelm->getAdrI();  // &VADI ????
      for(feInt i=0;i<NBRI;i++) {
        feInt I = (feInt) VADI[i];
        if(I<ordre) ddlNumberOfElements[I]++;
      }
    }
  }

  feInt totalNumberOfElements=0;
  for(feInt i=0;i<ordre;i++) totalNumberOfElements+=ddlNumberOfElements[i];
  // ================================================================
  // IDENTIFIER LES ELEMENTS ET LES FORMES BILINEAIES D'UN 
  // DEGRE DE LIBERTE -- PREPARARION DE L'ESPACE DE TRAVAIL
  // ================================================================
  ddlBiLinearForms  = new feInt* [ordre];
  ddlElms           = new feInt* [ordre];

 
  feInt* ptr = new feInt [totalNumberOfElements];
  for(feInt i=0;i<ordre;i++)     {
    ddlBiLinearForms[i]=ptr;
    ptr+=ddlNumberOfElements[i];
  }

  ptr = new feInt [totalNumberOfElements];
  for(feInt i=0;i<ordre;i++)    {
    ddlElms[i]=ptr;
    ptr+=ddlNumberOfElements[i];
  }

  // ================================================================
  // RECUPERER LES NUMEROS DES FORMES BIULINEAIRES ET DES ELEMENTS 
  // ================================================================
  for(feInt i=0;i<ordre;i++) ddlNumberOfElements[i]=0;
  for(feInt eq=0;eq<NumberOfBilinearForms;eq++)
  {    
    feBilinearForm* equelm = formMatrices[eq];
    feInt cncGeoTag = equelm->getCncGeoTag();
    feInt nbElems   = mesh->getNbElm(cncGeoTag);

    for(feInt el=0;el<nbElems;el++)
    {
      equelm->initialize_vadij_only(metaNumber, el);
      feInt             NBRI = equelm->getNiElm();
      std::vector<int>  VADI = equelm->getAdrI();

      for(feInt i=0;i<NBRI;i++) 
      {
        feInt I=VADI[i];
        if(I<ordre) 
        {
          feInt k=ddlNumberOfElements[I];
          ddlBiLinearForms   [I][k]=eq;
          ddlElms            [I][k]=el;
          ddlNumberOfElements[I]++;
        }
      }   
    }
  }
  // ================================================================
  // COMPTER LES COEFFICIENTS D'UNE RANGEE
  // ================================================================
  //feInt* liste        = new feInt [ordre];
  //bool*  ddl_rngcof   = new bool  [ordre];
  //feInt* ddl_nbrcof   = new feInt [ordre];
  liste               = new feInt [ordre];
  ddl_rngcof          = new bool  [ordre];
  feInt* ddl_nbrcof   = new feInt [ordre];
  feInt  eq    = 0;
  feInt  el    = 0;


  for(feInt i=0;i<ordre;i++)     {   
    ddl_rngcof[i]= false;
    liste[i]     = 0;
  }
  for(feInt I=0;I<ordre;I++)     {

    feInt nbrelm = ddlNumberOfElements[I];
    feInt nbrcof = 0;

  // ================================================================
  // Pour stoujours avoir un coefficient sur la diagonale
  // 18 fevrier 2011 mod -- AG
  // ================================================================ 
      ddl_rngcof[I]=true;
      liste[nbrcof]=I;  
      nbrcof++;

      for(feInt k=0;k<nbrelm;k++)     {  
        eq  = ddlBiLinearForms[I][k];
        el  = ddlElms         [I][k];

        feBilinearForm* equelm = formMatrices[eq];
        // feInt cncGeoTag = equelm->getCncGeoTag();
        equelm->initialize_vadij_only(metaNumber, el);
        feInt             NBRJ = equelm->getNiElm();
        std::vector<int>  VADJ = equelm->getAdrJ();

        for(feInt j=0;j<NBRJ;j++) {
          feInt J=VADJ[j];
          if(J<ordre) {
            if(!ddl_rngcof[J])  {
              ddl_rngcof[J]=true;              
              liste[nbrcof]=J;
              nbrcof++;
            }
          }
        }
      }

      ddl_nbrcof[I]=nbrcof;
      for(feInt i=0;i<nbrcof;i++) ddl_rngcof[liste[i]]=false;
  } 

  // ================================================================
  // ON CONSERVE NNZ
  // ================================================================

  for(feInt i=0;i<ordre;i++) nnz[i] = ddl_nbrcof[i];
  nz = 0;
  for(feInt i=0;i<ordre;i++) nz+=nnz[i];

  // ================================================================
  // LIBERER L'ESPACE DE TRAVAIL TEMPORAIRE
  // ================================================================
  //delete [] ddlNumberOfElements;

  // delete [] ddlBiLinearForms[0];
  // delete [] ddlBiLinearForms;
  // delete [] ddlElms[0];
  // delete [] ddlElms;

  // delete [] liste;
  // delete [] ddl_rngcof;
  delete [] ddl_nbrcof;
}

// ====================================================================
// Constructeur de la classe dérivée CSR MKL PARDISO
// ====================================================================

//========================================================================================================
// fonction pour comprarer deux entiers, pour le tri
// d'une liste d'entiers par la fonction qsort
//========================================================================================================
int compint (const void* a, const void* b ){
   const feInt* aa = (const feInt*) a;
   const feInt* bb = (const feInt*) b;
   int code=0;
   feInt val = (*aa-*bb);
   if      (val< 0) code = -1;
   else if (val==0) code =  0;
   else             code =  1;
   return code;
}

feCompressedRowStorageMklPardiso::feCompressedRowStorageMklPardiso(feMetaNumber *metaNumber, feMesh *mesh, std::vector<feBilinearForm*> &formMatrices)
: feCompressedRowStorage(metaNumber, mesh, formMatrices)
{

    // ================================================================
    // ALLOUER LA STRUCTURE CSR DE MKL PARDISO
    // ================================================================
    irangee = new feInt  [ordre];
    rangee  = new double [ordre];
    Ap      = new feInt  [ordre+1];
    Aj      = new feInt  [nz];   // EFMKLPARDISOint 
    //Ax      = new double[nz];
    // ================================================================
    // INITIALISATION
    // ================================================================
    for(feInt i=0;i<ordre  ;i++) rangee[i]=0;
    for(feInt i=0;i<ordre+1;i++) Ap[i]=0;
    for(feInt i=0;i<ordre  ;i++) Aj[i]=0;

    // =============================================================
    // GENERER LA STRUCTURE Ap 
    // =============================================================
    if(Ap!=NULL)
    {
        Ap[0]=0;
        for(feInt i=0;i<ordre;i++) 
           Ap[i+1]=Ap[i]+(feInt) nnz[i]; // EFMKLPARDISOint
    }  


    // ================================================================
    // COMPTER LES COEFFICIENTS D'UNE RANGEE
    // ON COMPTE DE NOUVEAU POUR CONSERVER LES NUMÉROS DES COLONNES
    // DANS LA STRUCTURE Aj
    // ================================================================
    feInt  eq           = 0;
    feInt  el           = 0;

    for(feInt i=0;i<ordre;i++)     {   
        ddl_rngcof[i]= false;
        liste[i]     = 0;
    }
    for(feInt I=0;I<ordre;I++)     {
        // Un pinteur sur le début de la rangée
        feInt* ptrAj = Aj + Ap[I];  

        feInt nbrelm = ddlNumberOfElements[I];
        feInt nbrcof = 0;

    // ================================================================
    // Pour stoujours avoir un coefficient sur la diagonale
    // 18 fevrier 2011 mod -- AG
    // ================================================================ 
        ddl_rngcof[I]=true;
        liste[nbrcof]=I;  
        nbrcof++;

        for(feInt k=0;k<nbrelm;k++)     {  
            eq  = ddlBiLinearForms[I][k];
            el  = ddlElms         [I][k];
    
            feBilinearForm* equelm = formMatrices[eq];
            // feInt cncGeoTag = equelm->getCncGeoTag();
            equelm->initialize_vadij_only(metaNumber, el);
            feInt             NBRJ = equelm->getNiElm();
            std::vector<int>  VADJ = equelm->getAdrJ();
 

            for(feInt j=0;j<NBRJ;j++) {
                feInt J=VADJ[j];
                if(J<ordre) {
                    if(!ddl_rngcof[J])  {
                        ddl_rngcof[J]=true;              
                        liste[nbrcof]=J;
                        nbrcof++;
                    }
                }
            }
        }

        for(feInt i=0;i<nbrcof;i++) ddl_rngcof[liste[i]]=false;

        // Il faut conserver la liste des colonnes
        qsort(liste, nbrcof, sizeof(feInt), compint);

        for(feInt i=0;i<nbrcof;i++)
        {
            ddl_rngcof[liste[i]]=false;
            *ptrAj=(feInt)liste[i];
            ptrAj++;
        }
    } 

    // ================================================================ 
    // PARDISO
    // ================================================================ 
    // ================================================================ 
    // CONVERSION BASE 0 C A BASE 1 FORTRAN
    // ================================================================ 
    for(feInt i=0;i<nz     ;i++) Aj[i]+=1;
    for(feInt i=0;i<ordre+1;i++) Ap[i]+=1;

}

void feCompressedRowStorageMklPardiso::print_info() {

    //for(feInt i=0;i<nz;i++) printf(" %ld \n", Aj[i]);
    for(feInt i=0;i<ordre;i++) {
         printf("rangee (%ld) ", i+1);
         for(feInt j=0;j<Ap[i+1]-Ap[i];j++) printf(" %ld ", Aj[Ap[i]-1+j]);
         printf("\n");   
        }

}

void  feCompressedRowStorageMklPardiso::matrixAddValues(double* Matrix, feInt nRow, feInt *Row, feInt nColumn, feInt *Column, double** Aij) {
    for(feInt i=0;i<nRow;i++) {
        feInt I = Row[i];
        if(I < ordre) {
            feInt debut   = Ap[I]-1;
            feInt fin     = Ap[I+1]-1;
            feInt ncf     = fin-debut;

            for (feInt j=0;j<ncf;j++) irangee[Aj[debut+j]-1]=debut+j;
            for (feInt j=0;j<nColumn;j++) 
            {
                feInt J = Column[j];
                if(J < ordre) 
                {   
                    printf("irangee %ld Aij %g\n",irangee[J], Aij[i][j]);
                    Matrix[irangee[J]] += Aij[i][j];
                }

            }
        }
    }
}



// void SOMMERSTRUCTUREP4          (feInt NBRI,const feInt* VADI,
//                                          feInt NBRJ,const feInt* VADJ,
//                                          const double** matriceelementaire,
//                      feInt ordre, 
//                      feInt* irangee,
//                                          feInt* trangee,
//                      feInt  numthread,
//                      EFMKLPARDISOint* Ap, EFMKLPARDISOint* Aj, double* Ax )
// {
//     for(feInt j=0;j<NBRJ;j++)
//     {
//         feInt J = VADJ[j]-1;
//         if( J < ordre ) trangee[J] = numthread;
//     }
//     for(feInt i=0;i<NBRI;i++)
//     {

//         feInt I = VADI[i]-1;
//         if( I < ordre )
//         {

//             // TRANSFERT DE LA MATRICE ELEMENTAIRE A UN VECTEUR RANGEE  DE DIMENSION
//             // EGALE A L'ORDRE DE LA MATRICE
//             // =====================================================================
//             // TRANSFERT DU VECTEUR RANGEE A LA STRUCTURE Ax
//             // =====================================================================
//             feInt debut   = Ap[I]-1;
//             feInt fin     = Ap[I+1]-1;
//             feInt ncf     = fin-debut;
//             //double* rangee= Ax+debut;
//             for (feInt j=0;j<ncf;j++)
//             {
//                 if(trangee[Aj[debut+j]-1]==numthread)
//                 irangee[Aj[debut+j]-1]=debut+j;
//             }
//             // METTRE A ZERO LES COEFFICIENTS DE LA COLONNE
//             // =====================================================================
//             for(feInt j=0;j<NBRJ;j++)
//             {
//                 feInt J = VADJ[j]-1;
//                 if(J < ordre) 
//                 {
//                     Ax[irangee[J]] += matriceelementaire[i][j];
//                 }
//             }
//         }
//     }
//     for(feInt j=0;j<NBRJ;j++)
//     {
//         feInt J = VADJ[j]-1;
//         if( J < ordre ) trangee[J] = -1;
//     }
// };





//void EF5_SYSEQS_MKLPARDISO::ALLOUERSTOCKAGEMATRICIEL(void)
//{
// ==========================================
// COMPTER LES ELEMENTS D'UN DEGRE DE LIBERTE
// ==========================================
//  feInt* ddl_nbrelm = new feInt [ordre];
//  for(feInt i=0;i<ordre;i++) ddl_nbrelm[i]=0;
//  feInt NombreEquation = prbequ->GETNBREQU();

//  for(feInt eq=0;eq<NombreEquation;eq++)
//  { 
//      EF5_EQUELM* equelm = prbequ->GETEQUELM(eq);
//      const EF5_PARTITIONGEOMETRIQUE* partition = equelm->GETPARTITION();

//      for(feInt el=0;el<partition->NOMBREDELEMENT();el++)
//      {
//              equelm->INITIALISERVAD(el);

//          feInt        NBRJ = equelm->GETNBREQUATION();
//          const feInt* VADJ = equelm->GETADRESSAGEEQUATION();
// // Est-ce qu'on devrait prendre les inconnues
//          for(feInt j=0;j<NBRJ;j++) 
//          {
//              feInt J=VADJ[j]-1;
//              if(J<ordre) ddl_nbrelm[J]++;
//          }
//      }
//  } 

//  feInt total=0;
//  for(feInt i=0;i<ordre;i++) total+=ddl_nbrelm[i];
// =============================================================
// IDENTIFIER LES ELEMENTS ET LES EQUATIONS D'UN DEGRE DE LIBERTE
// =============================================================
    // feInt** ddl_equation = new feInt* [ordre];
    // feInt** ddl_element  = new feInt* [ordre];
   
 //        //printf("EF5_SYSEQS ordre, total(%ld,%ld) \n",ordre,total); 
    // feInt* listeInit = new feInt [total];
    // feInt* ptr = listeInit;
    // for(feInt i=0;i<ordre;i++)
    // {
    //  ddl_equation[i]=ptr;
    //  ptr+=ddl_nbrelm[i];
    // }
    // listeInit = new feInt [total];
    // ptr = listeInit;
    // for(feInt i=0;i<ordre;i++)
    // {
    //  ddl_element[i]=ptr;
    //  ptr+=ddl_nbrelm[i];
    // }

//  for(feInt i=0;i<ordre;i++) ddl_nbrelm[i]=0;
//  for(feInt eq=0;eq<NombreEquation;eq++)
//  {    
// // ==========================================
// // COMPTER LES ELEMENTS D'UN DEGRE DE LIBERTE
// // ==========================================
//      EF5_EQUELM* equelm = prbequ->GETEQUELM(eq);
//      const EF5_PARTITIONGEOMETRIQUE* partition = equelm->GETPARTITION();

//      for(feInt el=0;el<partition->NOMBREDELEMENT();el++)
//      {
//              equelm->INITIALISERVAD(el);

//          feInt        NBRJ = equelm->GETNBREQUATION();
//          const feInt* VADJ = equelm->GETADRESSAGEEQUATION();
// // les inconnues?
//          for(feInt j=0;j<NBRJ;j++) 
//          {
//              feInt J=VADJ[j]-1;
//              if(J<ordre) 
//              {
//                  feInt k=ddl_nbrelm[J];
//                  ddl_equation[J][k]=eq;
//                  ddl_element [J][k]=el;
//                  ddl_nbrelm  [J]++;
//              }
//          }   
//      }
//  } 

// // ==========================================
// // COMPTER LES COEFFICIENTS D'UNE RANGEE
// // ==========================================
//  feInt* listeth      = new feInt [ordre*nbthreads];
//  feInt* ddl_rngcofth = new feInt [ordre*nbthreads];
//  feInt* ddl_nbrcof   = new feInt [ordre];
//  feInt numthread     = 0;
//  feInt eq            = 0;
//  feInt el            = 0;
//  feInt eqt           = 0;

// //DEBUT OPENMP
// #ifdef _OPENMP_  
// #pragma omp parallel shared(listeth,ddl_rngcofth)
// #endif
// {
// #ifdef _OPENMP_  
// #pragma omp for schedule(static)
// #endif
//  for(feInt i=0;i<ordre*nbthreads;i++) 
//  {   ddl_rngcofth[i]= FAUX;
//      listeth[i]     = 0;
//  }

// #ifdef _OPENMP_  
// #pragma omp for private(numthread,eq,el,eqt) schedule(dynamic)
// #endif
//  for(feInt J=0;J<ordre;J++) 
//  {
//      #ifdef _OPENMP_  
//      numthread = omp_get_thread_num();                   
//      #endif

//      feInt* liste = listeth + numthread*ordre;
//      feInt* ddl_rngcof = ddl_rngcofth + numthread*ordre;

//      feInt nbrelm = ddl_nbrelm[J];
//      feInt nbrcof = 0;


// // ==============================
// // pour s'assurer d'un coefficient sur la diagonale
// // 18 fevrier 2011 mod -- AG
// // ============================== 
// #ifdef DIAG 
//      ddl_rngcof[J]=VRAI;
//      liste[nbrcof]=J;  
//      nbrcof++;
// #endif


//      for(feInt k=0;k<nbrelm;k++)
//      {
//          eq  = ddl_equation[J][k];
//          el  = ddl_element [J][k];
//          eqt = eq +numthread*NombreEquation;

//          EF5_EQUELM* equelm = prbequ->GETEQUELM(eqt);

//          equelm->INITIALISERVAD(el);
 
//          feInt        NBRJ = equelm->GETNBRINCONNUE();
//          const feInt* VADJ = equelm->GETADRESSAGEINCONNUE();

//          for(feInt i=0;i<NBRJ;i++)
//          {
//              feInt I=VADJ[i]-1;
//              if(I<ordre) {
//                  if(!ddl_rngcof[I])
//                  {
//                      ddl_rngcof[I]=VRAI;              
//                      liste[nbrcof]=I;
//                      nbrcof++;
//                  }
//              }
//          }
//      }

//  ========= ICI ====== pour la conversion

//      ddl_nbrcof[J]=nbrcof;
//      for(feInt i=0;i<nbrcof;i++)
//          ddl_rngcof[liste[i]]=FAUX;
//  } 
    
// #ifdef _OPENMP_  
// #pragma omp single
// #endif
//     {
//  feInt somme =0;
//  for(feInt i=0;i<ordre;i++) somme += ddl_nbrcof[i];

//  nz =  (EFMKLPARDISOint) somme;
// // ===========================
// // ALLOUER L'ESPACE DE TRAVAIL
// // ===========================
//  rangee  = new double [ordre];
//  irangee = new feInt  [ordre];
//         trangee = new feInt  [ordre];
//  Ap      = new EFMKLPARDISOint [ordre+1];
//  Aj      = new EFMKLPARDISOint [nz];
//  Ax      = new double[nz];
// // ==============
// // INITIALISATION
// // ==============
//  for(feInt i=0;i<ordre  ;i++) rangee[i]=0;
//  for(feInt i=0;i<ordre  ;i++) irangee[i]=-1;
//  for(feInt i=0;i<ordre  ;i++) trangee[i]=-1;
//  for(feInt i=0;i<ordre+1;i++) Ap[i]=0;
// // =============================================================
// // GENERER LA STRUCTURE Ap 
// // =============================================================
//  if(Ap!=NULL)
//  {
//      Ap[0]=0;
//      for(feInt i=0;i<ordre;i++) 
//            Ap[i+1]=Ap[i]+(EFMKLPARDISOint) ddl_nbrcof[i];
//  }
//     }

// #ifdef _OPENMP_  
// #pragma omp for schedule(static)
// #endif
//  for(feInt i=0;i<nz   ;i++) 
//  {
//      Aj[i]= 0; 
//      Ax[i]= 0;
//  }

// #ifdef _OPENMP_  
// #pragma 




//  for(feInt J=0;J<ordre;J++)
//  {
//         #ifdef _OPENMP_  
//      numthread = omp_get_thread_num();                   
//      #endif

//         feInt* liste = listeth + numthread*ordre;
//      feInt* ddl_rngcof = ddl_rngcofth + numthread*ordre;

//      EFMKLPARDISOint* ptrAj = Aj + Ap[J];
//      feInt nbrelm = ddl_nbrelm[J];
//      EFMKLPARDISOint nbrcof = 0;
// // ==============================
// // pour s'assurer d'un coefficient sur la diagonale
// // 18 fevrier 2011 mod -- AG
// // ============================== 
// #ifdef DIAG
//      ddl_rngcof[J]=VRAI;
//      liste[nbrcof]=J;
//      nbrcof++;
// #endif       

//      for(feInt k=0;k<nbrelm;k++)
//      {
//          feInt eq=ddl_equation[J][k];
//          feInt el=ddl_element [J][k];
//          eqt = eq +numthread*NombreEquation;

//          EF5_EQUELM* equelm = prbequ->GETEQUELM(eqt);

//          equelm->INITIALISERVAD(el);
 
//          feInt        NBRJ = equelm->GETNBRINCONNUE();
//          const feInt* VADJ = equelm->GETADRESSAGEINCONNUE();

//          for(feInt i=0;i<NBRJ;i++)
//          {
//              feInt I=VADJ[i]-1;
//              if(I<ordre) {
//                  if(!ddl_rngcof[I])
//                  {
//                      ddl_rngcof[I]=VRAI;
//                      liste[nbrcof]=I;
//                      nbrcof++;
//                  }
//              }
//          }

//      }
//      qsort(liste, nbrcof, sizeof(feInt), compint);

//      for(feInt i=0;i<nbrcof;i++)
//      {
//          ddl_rngcof[liste[i]]=FAUX;
//          *ptrAj=(EFMKLPARDISOint)liste[i];
//          ptrAj++;
//      }
//  } 



// // ===========================
// // PARDISO
// // =========================== 
// // ===========================
// // CONVERSION BASE 0 C A BASE 1 FORTRAN
// // ===========================
// #ifdef _OPENMP_  
// #pragma omp for schedule(static)
// #endif
//  for(feInt i=0;i<nz     ;i++) Aj[i]+=1;

// } //Fin OPENMP
//  for(feInt i=0;i<ordre+1;i++) Ap[i]+=1;

// // ===========================
// // LIBERER L'ESPACE DE TRAVAIL
// // ===========================
//  delete [] listeth;
//  delete [] ddl_rngcofth;
//  delete [] ddl_nbrelm;
//  if(ordre>0) delete [] ddl_equation[0];
//  delete [] ddl_equation;
//  if(ordre>0) delete [] ddl_element[0];
//  delete [] ddl_element;

// // ============================
// // TEST
// // ============================
//  //printf("ordre = %ld, nz = %ld \n",ordre, nz);
//  //exit(-1);
// };
