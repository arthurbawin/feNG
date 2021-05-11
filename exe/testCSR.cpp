#include "feNG.h"
#include "feFunction.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feCncGeo.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feBilinearForm.h"
#include "feLinearSystem.h"
#include "feSolver.h"
#include "feCompressedRowStorage.h"

double f1(const double t, const std::vector<double> x, const std::vector<double> par){
	return x[0];
}

double f2(const double t, const std::vector<double> x, const std::vector<double> par){
	return 0.0;
}


int main( int argc, char** argv )
{
	int ret = 0;

	// ============================================================
	// FONCTIONS
	// ============================================================
	 feFunction* feSolex = new feFunction(f1,{0,1});

	// ============================================================
	// MAILLAGE
	// ============================================================
	double xa    = 0;
	double xb    = 5;
	int    nbelm = 3;
    feMesh1DP1 mesh(xa,xb,nbelm,"BXA","BXB","M1D");
    //femaillage.printInfo();
	// ============================================================
	// FESPACE
	// ============================================================
  	feSpace1DP0 U_BXA (&mesh, "U", "BXA", feSolex);
  	feSpace1DP0 U_BXB (&mesh, "U", "BXB", feSolex);
  	feSpace1DP2 U_M1D (&mesh, "U", "M1D", feSolex);
  	feSpace1DP0 V_BXA (&mesh, "V", "BXA", feSolex);
  	feSpace1DP0 V_BXB (&mesh, "V", "BXB", feSolex);
  	feSpace1DP2 V_M1D (&mesh, "V", "M1D", feSolex);
  	// feSpace1DP0 V_BXA = feSpace1DP0(&mesh, "V", "BXA", feSolex);
  	// feSpace1DP0 V_BXB = feSpace1DP0(&mesh, "V", "BXB", feSolex);
  	// feSpace1DP3 V_M1D = feSpace1DP3(&mesh, "V", "M1D", feSolex);
  	std::vector<feSpace*> fespace = {&U_BXA, &U_BXB, &U_M1D,&V_BXA, &V_BXB, &V_M1D};
  	std::vector<feSpace*> feEssBC = {&U_BXA, &U_BXB,&V_BXA, &V_BXB};
	// ============================================================
	// NUMEROTATION
	// ============================================================
	feMetaNumber *metaNumber = new feMetaNumber(&mesh, fespace, feEssBC);

	// ============================================================
	// FORMES BILINEAIRES
	// ============================================================
	feSysElm_1D_Diffusion *diffusion1D = new feSysElm_1D_Diffusion(1.0, nullptr);
  	std::vector<feSpace*> spaceDiffusion1D_U_M1D = {&U_M1D};
  	std::vector<feSpace*> spaceDiffusion1D_V_M1D = {&V_M1D};

  	int nQuadraturePoints = 3; // TODO : change to deg
  	feBilinearForm *bDiff_U_M1D = new feBilinearForm(spaceDiffusion1D_U_M1D, &mesh, nQuadraturePoints, diffusion1D);
  	feBilinearForm *bDiff_V_M1D = new feBilinearForm(spaceDiffusion1D_V_M1D, &mesh, nQuadraturePoints, diffusion1D);
  	//feBilinearForm *bDiff_U_M1D = new feBilinearForm({&U_M1D}, &mesh, nQuadraturePoints, diffusion1D);

  	std::vector<feBilinearForm*> formMatrices  = {bDiff_U_M1D,bDiff_V_M1D};


	//feCompressedRowStorage feCRS(metaNumber, &mesh, fespace);
	feCompressedRowStorage* feCRS = new feCompressedRowStorage(metaNumber, &mesh, formMatrices);
	int    nz = feCRS->getNz();
	feInt* nnz= feCRS->getNnz();
	printf("nz = %d \n",nz);
	for(int i=0;i<metaNumber->getNbUnknowns();i++) {
		printf("rangée (%d) nz(%ld) \n", i, nnz[i]);

	}
	// ============================================================
	// SOLUTION
	// ============================================================
    delete feSolex;
    delete metaNumber;
   	return ret;
}
