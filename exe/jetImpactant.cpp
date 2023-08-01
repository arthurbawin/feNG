
#include "feAPI.h"

static double geo_D = 1.;
static double geo_h = 20.;
static double geo_H =  4. * geo_D;
static double geo_L = 15. * geo_D;

bool inInjectionPipe(double* x, bool strictly = false)
{
  double xmin = geo_L - geo_D/2.;
  double xmax = geo_L + geo_D/2.;
  double ymin = geo_H;
  double ymax = geo_h;
  if(strictly) {
    // Treat the inside at the connection of the two pipes
    return (x[0] > xmin) && (x[0] < xmax) && (x[1] >= ymin) && (x[1] < ymax);
  } else
    return (x[0] >= xmin) && (x[0] <= xmax) && (x[1] >= ymin) && (x[1] <= ymax);
}

bool inBottomRectangle(double* x, bool strictly = false)
{
  if(strictly) {
    // Treat the inside at the connection of the two pipes
    if(geo_L - geo_D/2. < x[0] && x[0] < geo_L + geo_D/2.)
      return (x[0] > 0.) && (x[0] < 2.*geo_L) && (x[1] > 0.) && (x[1] <= geo_H);
    else
      return (x[0] > 0.) && (x[0] < 2.*geo_L) && (x[1] > 0.) && (x[1] <  geo_H);
  } else
    return (x[0] >= 0.) && (x[0] <= 2.*geo_L) && (x[1] >= 0.) && (x[1] <= geo_H);
}

bool inside(double* x, bool strictly = false)
{
  return inInjectionPipe(x, strictly) || inBottomRectangle(x, strictly);
}

void vInlet(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  double x = pos[0];
  double l1 = par[0];
  double D = par[1];
  double Re = par[2];

  // double uMax = (D*D - l1*l1)/4.;

  res[0] = 0.;
  res[1] = - Re/2. * (x-l1) * (D - (x-l1));
  if(res[1] > 0) {
    feInfo("y-velocity is %f", res[1]);
    exit(-1);
  }
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);

  int orderVelocity = 2;
  int orderPressure = 1;
  int degreeQuadrature = 8;
  bool curveTheMesh = false;

  setVerbose(1);

  feConstantFunction zero(0.);
  feConstantFunction one(1.);
  feConstantFunction minusOne(-1.);
  feConstantVectorFunction zeroVector({0., 0.});

  #define NMESH 1
  #define RESTART_AT 0
  int nMesh = NMESH;
  int TARGETN[7] = {500, 1000, 2000, 4000, 8000, 16000, 32000};
  int targetN[NMESH], nV[NMESH];
  double E[NMESH];

  for(int i = 0; i < nMesh; ++i) {
    targetN[i] = TARGETN[ RESTART_AT + i ];
  }

  std::string meshFile;

  feSolutionContainer containerN;

  for(int iV = 0; iV < nMesh; ++iV)
  {
    feSolutionContainer *container = nullptr;

    double theT0 = 0.;
    double dt = 0.05;
    int nTimeSteps = 10;
    int currentStep = 0;
    int nAdaptationCycles = 100;

    double Re = 50.;

    double D = 1.;
    double l1 = 15. * D - D/2.;
    feVectorFunction inlet(vInlet, {l1, D, Re});


    bool adapt = true;
    bool curved = true;

    std::string outputDirectory = "jetImpactant/";

    for(int iAdapt = 0; iAdapt < nAdaptationCycles; ++iAdapt)
    {
      if(adapt && iAdapt > 0) {
        meshFile = "adapted.msh";
      } else {
        if(curved) {
          meshFile = "../data/jetImpactant_P2.msh";
          // meshFile = "ref_adapted.msh";
          // meshFile = "ref_curved.msh";
	      } else {
          meshFile = "../data/jetImpactant.msh";
        }
      }

      // ===== Solver starts here ==========================
      feMesh2DP1 mesh(meshFile, curved);
      nV[iV] = mesh.getVertices().size();

      feConstantFunction density(1.);
      feConstantFunction viscosity(D/Re);

      feSpace *uDomaine = nullptr, *pDomaine = nullptr, *uIn = nullptr, *uNoSlip = nullptr, *uBottom = nullptr, *wDomaine = nullptr;
      feCheck(createFiniteElementSpace(           uIn, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",     "Inlet", degreeQuadrature, &inlet));
      feCheck(createFiniteElementSpace(       uBottom, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_Tw", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace(       uNoSlip, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_T0", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace(      uDomaine, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",   "Domaine", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace(      pDomaine, &mesh, elementType::LAGRANGE,        orderPressure, "P",   "Domaine", degreeQuadrature, &zero));
      feCheck(createFiniteElementSpace(wDomaine, &mesh, elementType::LAGRANGE,        orderVelocity, "W",  "Domaine", degreeQuadrature, &zero));
      std::vector<feSpace*> spaces = {uIn, uBottom, uNoSlip, uDomaine, pDomaine, wDomaine};
      std::vector<feSpace*> essentialSpaces = {uNoSlip, uIn, uBottom};

      feMetaNumber numbering(&mesh, spaces, essentialSpaces);
      feSolution sol(numbering.getNbDOFs(), spaces, essentialSpaces);

      if(iAdapt > 0) {
        sol.setSolFromContainer(container, 0);
        for(auto *space : spaces)
          space->setDOFInitialization(dofInitialization::PREVIOUS_SOL);
        uIn->setDOFInitialization(dofInitialization::NODEWISE);
      }

      feBilinearForm *diffU = nullptr, *gradP = nullptr, *divU = nullptr, *conv = nullptr, *supg = nullptr, *trmass = nullptr;
      feBilinearForm *vorticity = nullptr, *massvorticity = nullptr;
      feCheck(createBilinearForm( gradP, {uDomaine, pDomaine}, new feSysElm_MixedGradient(&minusOne)               ));
      feCheck(createBilinearForm( diffU,           {uDomaine}, new feSysElm_VectorDiffusion(&one, &viscosity)      ));
      feCheck(createBilinearForm(  divU, {pDomaine, uDomaine}, new feSysElm_MixedDivergence(&one)                  ));
      feCheck(createBilinearForm(  conv,           {uDomaine}, new feSysElm_VectorConvectiveAcceleration(&minusOne)));
      feCheck(createBilinearForm(trmass,           {uDomaine}, new feSysElm_TransientVectorMass(-1.)               ));
      // feCheck(createBilinearForm(  supg, {uDomaine, pDomaine}, new feSysElm_NS_SUPG_PSPG(&one, &density, &viscosity, &zeroVector) ));
      
      // Vorticity equation
      feCheck(createBilinearForm( vorticity, {wDomaine, uDomaine}, new feSysElm_MixedCurl(&minusOne) ));
      feCheck(createBilinearForm( massvorticity, {wDomaine}, new feSysElm_Mass(&one) ));

      std::vector<feBilinearForm*> forms = {gradP, diffU, divU, conv, trmass};

      if(wDomaine) {
        forms.push_back(vorticity);
        forms.push_back(massvorticity);
      }

      feLinearSystem *linearSystem;
      feCheck(createLinearSystem(linearSystem, MKLPARDISO, forms, numbering.getNbUnknowns(), argc, argv));

      // curveTheMesh = iAdapt > 0;
      curveTheMesh = 0;

      feExporter *exporter;
      feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
      int exportEveryNSteps = 1;
      std::string vtkFileRoot = outputDirectory + "jetImpactant_nv" + std::to_string(targetN[iV]) + "_";
      if(!adapt)
        vtkFileRoot = vtkFileRoot + "noAdapt_";
      if(adapt && curveTheMesh)
        vtkFileRoot = vtkFileRoot + "curved_";
      feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

      feNorm *norm;
      feCheck(createNorm(norm, VECTOR_L2, {uDomaine}, &sol));
      std::vector<feNorm*> norms = {norm};

      feTolerances tol{1e-10, 1e-10, 200};

      TimeIntegrator *solver;

      double t0 = theT0;
      double t1 = t0 + dt*nTimeSteps;
      feCheck(createTimeIntegrator(solver, BDF2, tol, linearSystem, &numbering, &sol, &mesh, norms, exportData, t0, t1, nTimeSteps));

      if(iAdapt > 0) {
        solver->setCurrentStep(currentStep);
      }

      // Solve
      for(int ii = 0; ii < nTimeSteps; ++ii) {
        feCheck(solver->makeStep());
      }

      currentStep = solver->getCurrentStep();
      theT0 = sol.getCurrentTime();
      container = solver->getSolutionContainer();

      // ===== Solver ends here ==========================

      if(adapt) {
        // ===== Mesh adaptation starts here ==========================
        feMetricOptions options(meshFile);
        options.nTargetVertices = targetN[iV];
        options.LpNorm = 100.;
        options.insideCallback = inside;
        options.hMin = 1e-10;
        options.hMax = 10.;
          
        if(orderVelocity == 1) {
          options.method = adaptationMethod::ANISO_P1;
          options.enableGradation = true;
          options.gradation = 1.5;
        } else {
          // if(curveTheMesh) {
          //   options.method = adaptationMethod::CURVED_EXTREME_SIZES;
          //   options.enableGradation = true;
          //   options.gradation = 1.3;
          //   // options.eTargetError = 1e-3;
          // } else {
            options.method = adaptationMethod::ANISO_PN;
            options.enableGradation = true;
            options.gradation = 1.3;
            options.logSimplexOptions.nThetaPerQuadrant = 20;
            options.logSimplexOptions.tol = 0.02;
            options.logSimplexOptions.uniformErrorCurve = true;
            options.logSimplexOptions.numLoopsUniformErrorCurve = 5;
          // }
        }

        bool isBackMeshP2 = true;
        bool setGmshModelToP1 = true; // Set this to false only to test convergence of metric interpolation
        bool curveMMGmesh = true;
        bool reconstructAtHighOrderNodes = curved;

        // ============================================================================================================================================================
        // Step 1: Generate aniso mesh with MMG
        // feNewRecovery recV(uDomaine, 1, &mesh, &sol, meshFile, "recoveredDerivativesV.msh", reconstructAtHighOrderNodes, false, nullptr, &numbering);
        // std::vector<feNewRecovery*> recoveredFields = {&recV};
        feNewRecovery recW(wDomaine, 0, &mesh, &sol, meshFile, "recoveredDerivativesW.msh", reconstructAtHighOrderNodes, false, nullptr, &numbering);
        std::vector<feNewRecovery*> recoveredFields = {&recW};
        feCheck(mesh.adapt(recoveredFields, options, spaces, essentialSpaces, uDomaine,
          nullptr, nullptr, nullptr, curveTheMesh, isBackMeshP2, setGmshModelToP1, curveMMGmesh, curveToMinimize::LENGTH));

        // ============================================================================================================================================================
        // Interpolate solution on adapted mesh
        feMesh2DP1 nextMesh("adapted.msh", curved);
        feSpace *uDomaine = nullptr, *pDomaine = nullptr, *uIn = nullptr, *uNoSlip = nullptr, *uBottom = nullptr, *wDomaine = nullptr;
        feCheck(createFiniteElementSpace(     uIn, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",     "Inlet", degreeQuadrature, &inlet));
        feCheck(createFiniteElementSpace( uBottom, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_Tw", degreeQuadrature, &zeroVector));
        feCheck(createFiniteElementSpace( uNoSlip, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_T0", degreeQuadrature, &zeroVector));
        feCheck(createFiniteElementSpace(uDomaine, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",   "Domaine", degreeQuadrature, &zeroVector));
        feCheck(createFiniteElementSpace(pDomaine, &nextMesh, elementType::LAGRANGE,        orderPressure, "P",   "Domaine", degreeQuadrature, &zero));
        feCheck(createFiniteElementSpace(wDomaine, &nextMesh, elementType::LAGRANGE,        orderVelocity, "W",   "Domaine", degreeQuadrature, &zero));
        std::vector<feSpace*> nextSpaces = {uIn, uBottom, uNoSlip, uDomaine, pDomaine, wDomaine};
        std::vector<feSpace*> nextEssentialSpaces = {uNoSlip, uIn, uBottom};
        feMetaNumber nextNumbering(&nextMesh, nextSpaces, nextEssentialSpaces);
        feSolution nextSol(nextNumbering.getNbDOFs(), nextSpaces, nextEssentialSpaces);
        feCheck(mesh.transfer(&nextMesh, &numbering, &nextNumbering, container, spaces, essentialSpaces, nextSpaces));
        nextSol.setSolFromContainer(container, 0);
        feExporter *nextExporter;
        feCheck(createVisualizationExporter(nextExporter, VTK, &nextNumbering, &nextSol, &nextMesh, nextSpaces));
        feCheck(nextExporter->writeStep("afterTransfer.vtk"));

        // ============================================================================================================================================================
        // Step 2: Curve the MMG mesh
        // curveMMGmesh = true;
        // feNewRecovery recW2(wDomaine, 0, &nextMesh, &nextSol, "adapted.msh", "recoveredDerivativesW.msh", reconstructAtHighOrderNodes, false, nullptr, &nextNumbering);
        // recoveredFields[0] = &recW2;

        // curveTheMesh = 1;
        // options.backgroundMeshfile = "adapted.msh";
        // feCheck(mesh.adapt(recoveredFields, options, nextSpaces, nextEssentialSpaces, uDomaine,
        //   nullptr, nullptr, nullptr, curveTheMesh, isBackMeshP2, setGmshModelToP1, curveMMGmesh, curveToMinimize::LENGTH));

        // ============================================================================================================================================================
        // Interpolate solution on curved mesh
        // feMesh2DP1 curvedMesh("adapted.msh", curved);
        // feSpace *uDomaine2 = nullptr, *pDomaine2 = nullptr, *uIn2 = nullptr, *uNoSlip2 = nullptr, *uBottom2 = nullptr, *wDomaine2 = nullptr;
        // feCheck(createFiniteElementSpace(     uIn2, &curvedMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",     "Inlet", degreeQuadrature, &inlet));
        // feCheck(createFiniteElementSpace( uBottom2, &curvedMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_Tw", degreeQuadrature, &zeroVector));
        // feCheck(createFiniteElementSpace( uNoSlip2, &curvedMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "NoSlip_T0", degreeQuadrature, &zeroVector));
        // feCheck(createFiniteElementSpace(uDomaine2, &curvedMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",   "Domaine", degreeQuadrature, &zeroVector));
        // feCheck(createFiniteElementSpace(pDomaine2, &curvedMesh, elementType::LAGRANGE,        orderPressure, "P",   "Domaine", degreeQuadrature, &zero));
        // feCheck(createFiniteElementSpace(wDomaine2, &curvedMesh, elementType::LAGRANGE,        orderVelocity, "W",   "Domaine", degreeQuadrature, &zero));
        // std::vector<feSpace*> curvedSpaces = {uIn2, uBottom2, uNoSlip2, uDomaine2, pDomaine2, wDomaine2};
        // std::vector<feSpace*> curvedEssentialSpaces = {uNoSlip2, uIn2, uBottom2};
        // feMetaNumber curvedNumbering(&curvedMesh, curvedSpaces, curvedEssentialSpaces);
        // feSolution curvedSol(curvedNumbering.getNbDOFs(), curvedSpaces, curvedEssentialSpaces);
        // feCheck(nextMesh.transfer(&curvedMesh, &nextNumbering, &curvedNumbering, container, nextSpaces, nextEssentialSpaces, curvedSpaces));
        // curvedSol.setSolFromContainer(container, 0);
        // feExporter *curvedExporter;
        // feCheck(createVisualizationExporter(curvedExporter, VTK, &curvedNumbering, &curvedSol, &curvedMesh, curvedSpaces));
        // feCheck(curvedExporter->writeStep("afterTransferOnCurvedMesh.vtk"));

        // ===== Mesh adaptation and error estimation end here =================================
      } // if adapt
    } // for iAdapt
  } // for iV

  petscFinalize();
  return 0;
}
