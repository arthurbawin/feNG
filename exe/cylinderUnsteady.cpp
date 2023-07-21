
#include "feAPI.h"

double BBOX_SIZE = 1.;
static double xmin = -BBOX_SIZE;
static double xmax =  BBOX_SIZE;
static double ymin = -BBOX_SIZE;
static double ymax =  BBOX_SIZE;

bool inside(double* x, bool strictly = false){
  feInfo("Modifiy inside for this geometry");
  exit(-1);
  if(strictly){
    return (x[0] > xmin) && (x[1] > ymin) && (x[0] < xmax) && (x[1] < ymax);
  } else{
    return (x[0] >= xmin) && (x[1] >= ymin) && (x[0] <= xmax) && (x[1] <= ymax);
  }
}

void u1(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  double y = pos[1];
  double h = par[0];
  double t0 = par[1];

  // double epsilon = (y > h/2.)  ? (1. - fmin(t/t0, 1.)) : 0.;
  // double epsilon = (y >= h/2.) ? sin(2.*M_PI*y/h) : 0.;
  double epsilon = (y >= h/2.) ? 0.25 : 0.;

  // res[0] = 1. + (1. * y/h) * (1. - fmin(t/t0, 1.));
  // res[0] = 1. + fabs(epsilon) * (1. - fmin(t/t0, 1.));
  // res[0] = 1. + fabs(epsilon);
  res[0] = 0.25 + y/h;
  // res[0] = 1. + (1. * y/h) * fmax(0., 1. - t/t0);
  res[1] = 0.;
}

void u2(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  res[0] = 1.;
  res[1] = 0.;
}

void solU(const double t, const std::vector<double> &pos, const std::vector<double> &par, std::vector<double> &res)
{
  double y = pos[1];
  double H = par[0];
  double dpdx = par[1];

  res[0] = -dpdx/2. * y * (H - y);
  res[1] = 0.;
}

int main(int argc, char **argv)
{
  petscInitialize(argc, argv);

  int orderVelocity = 2;
  int orderPressure = 1;
  int degreeQuadrature = 15;
  bool curveTheMesh = false;
  bool appendToErrorFile = false;

  setVerbose(1);

  feConstantFunction zero(0.);
  feConstantFunction one(1.);
  feConstantFunction minusOne(-1.);
  feConstantVectorFunction zeroVector({0., 0.});

  #define NMESH 1
  #define RESTART_AT 2
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

    feVectorFunction inletUniform(u2);

    // Start the adaptation loop
    // int nAdaptationCycles = 200 / (nTimeSteps * dt);
    int nAdaptationCycles = 1;
    int nStepsToIncreaseRe = nAdaptationCycles / 2;

    double D = 1.;
    double h = 16.;
    double T0 = 2.; // Time integration starts at 0 from steady solution
    feVectorFunction inlet(u1, {h, T0});
    double dpdx = -1.0;
    feVectorFunction uExact(solU, {h, dpdx});

    double Re = 1.;
    double Re_final = 200.;
    double dRe = (Re_final - Re) / nStepsToIncreaseRe;

    bool adapt = false;
    bool curved = false;

    std::vector<double> forces_x(nAdaptationCycles);
    std::vector<double> forces_y(nAdaptationCycles);

    for(int iAdapt = 0; iAdapt < nAdaptationCycles; ++iAdapt)
    {
      if(adapt && iAdapt > 0) {
        meshFile = "anisoadapted.msh";
      } else {
        if(curved) {
          meshFile = "../data/cylindre_P2.msh";
        } else {
          // meshFile = "../data/cylindre_P1.msh";
          meshFile = "../data/poiseuille.msh";
        }
      }

      // meshFile = "../data/cylindreStructured.msh";
      // meshFile = "../data/VonKarmanV4.msh";

      // meshFile = "../data/cylindreP2_" + std::to_string(iAdapt+1) + ".msh";

      // ===== Solver starts here ==========================
      feMesh2DP1 mesh(meshFile, curved);
      nV[iV] = mesh.getVertices().size();

      feConstantFunction density(1.);
      feConstantFunction viscosity(D/Re);

      feSpace *uDomaine = nullptr, *pDomaine = nullptr, *uIn = nullptr, *uTop = nullptr, *uBottom = nullptr, *uNoSlip = nullptr, *wDomaine = nullptr;
      feSpace *pCylindre = nullptr;
      feSpace *lambdaCylindre = nullptr;
      feSpace *uDroite = nullptr;

      // Cylindre
      // feCheck(createFiniteElementSpace( uIn, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",   "Gauche", degreeQuadrature, &inletUniform));
      // feCheck(createFiniteElementSpace(    uTop, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",     "Haut", degreeQuadrature, &zeroVector));
      // feCheck(createFiniteElementSpace( uBottom, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",      "Bas", degreeQuadrature, &zeroVector));
      // feCheck(createFiniteElementSpace( uNoSlip, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Cylindre", degreeQuadrature, &zeroVector));
      // feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Domaine", degreeQuadrature, &inlet));
      // feCheck(createFiniteElementSpace(pDomaine, &mesh, elementType::LAGRANGE,        orderPressure, "P",  "Domaine", degreeQuadrature, &zero));
      // feCheck(createFiniteElementSpace(pCylindre, &mesh, elementType::LAGRANGE,        orderPressure, "P",  "Cylindre", degreeQuadrature, &zero));
      // // feCheck(createFiniteElementSpace(wDomaine, &mesh, elementType::LAGRANGE,        orderVelocity, "W",  "Domaine", degreeQuadrature, &zero));
      // // feCheck(createFiniteElementSpace(lambdaCylindre, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "L",  "Cylindre", degreeQuadrature, &zeroVector));
      // std::vector<feSpace*> spaces = {uTop, uBottom, uNoSlip, uDomaine, pDomaine, pCylindre, uIn};
      // std::vector<feSpace*> essentialSpaces = {uIn, uNoSlip};
      // uTop->setEssentialComponent(1, true);
      // uBottom->setEssentialComponent(1, true);

      // Poiseuille
      feCheck(createFiniteElementSpace(     uIn, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Gauche", degreeQuadrature, &inletUniform));
      feCheck(createFiniteElementSpace( uDroite, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Droite", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace( uNoSlip, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "NoSlip", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace(uDomaine, &mesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domaine", degreeQuadrature, &zeroVector));
      feCheck(createFiniteElementSpace(pDomaine, &mesh, elementType::LAGRANGE,        orderPressure, "P", "Domaine", degreeQuadrature, &zero));
      std::vector<feSpace*> spaces = {uIn, uDroite, uNoSlip, uDomaine, pDomaine};
      std::vector<feSpace*> essentialSpaces = {uIn, uNoSlip};
      // uDroite->setEssentialComponent(1, true);   
      // uIn->setEssentialComponent(0, true);   

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
      feBilinearForm *massUL = nullptr, *massLU = nullptr, *srcLambda = nullptr;
      feBilinearForm *divSigma = nullptr;
      feCheck(createBilinearForm( gradP, {uDomaine, pDomaine}, new feSysElm_MixedGradient(&minusOne)                    ));
      feCheck(createBilinearForm( diffU,           {uDomaine}, new feSysElm_VectorDiffusion(&one, &viscosity) ));
      feCheck(createBilinearForm( divSigma, {uDomaine, pDomaine}, new feSysElm_DivergenceNewtonianStress(&minusOne, &viscosity) ));
      feCheck(createBilinearForm(  divU, {pDomaine, uDomaine}, new feSysElm_MixedDivergence(&one)                  ));
      feCheck(createBilinearForm(  conv,           {uDomaine}, new feSysElm_VectorConvectiveAcceleration(&minusOne)));
      feCheck(createBilinearForm(trmass,           {uDomaine}, new feSysElm_TransientVectorMass(-1.)));
      // feCheck(createBilinearForm(  supg, {uDomaine, pDomaine}, new feSysElm_NS_SUPG_PSPG(&one, &density, &viscosity, &zeroVector) ));
      
      // Vorticity equation
      // feCheck(createBilinearForm( vorticity, {wDomaine, uDomaine}, new feSysElm_MixedCurl(&minusOne) ));
      // feCheck(createBilinearForm( massvorticity, {wDomaine}, new feSysElm_Mass(&one) ));

      // Enforce no-slip on the cylinder with a Lagrange multiplier
      // feConstantVectorFunction cylinderVelocity({0.,0.});
      // feCheck(createBilinearForm(   massUL, {uNoSlip, lambdaCylindre}, new feSysElm_MixedVectorMass(&one) ));
      // feCheck(createBilinearForm(   massLU, {lambdaCylindre, uNoSlip}, new feSysElm_MixedVectorMass(&one) ));
      // feCheck(createBilinearForm(srcLambda,          {lambdaCylindre}, new feSysElm_VectorSource(&cylinderVelocity)));

      // std::vector<feBilinearForm*> forms = {gradP, diffU, divU, conv, trmass, massUL, massLU, srcLambda};
      // std::vector<feBilinearForm*> forms = {gradP, diffU, divU, conv, trmass};
      std::vector<feBilinearForm*> forms = {gradP, diffU, divU, conv, trmass};
      // std::vector<feBilinearForm*> forms = {divSigma, divU, conv, trmass};

      feLinearSystem *linearSystem;
      // feCheck(createLinearSystem(linearSystem, PETSC, forms, numbering.getNbUnknowns(), argc, argv));
      feCheck(createLinearSystem(linearSystem, MKLPARDISO, forms, numbering.getNbUnknowns()));
      // feCheck(createLinearSystem(linearSystem, PETSC, forms, numbering.getNbDOFs(), argc, argv));
      // linearSystem->setDisplayRHSInConsole(true);

      feExporter *exporter;
      feCheck(createVisualizationExporter(exporter, VTK, &numbering, &sol, &mesh, spaces));
      int exportEveryNSteps = 1;
      // std::string vtkFileRoot = "cylinder_nv_" + std::to_string(targetN[iV]) + "_iter_" + std::to_string(iAdapt) + "_Re_" + std::to_string((int) floor(Re)) + "_";
      std::string vtkFileRoot = "cylinder_nv_" + std::to_string(targetN[iV]) + "_";
      if(!adapt)
        vtkFileRoot = vtkFileRoot + "noAdapt_";
      feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};

      sol.initializeUnknowns(&mesh);
      sol.initializeEssentialBC(&mesh);
      feCheck(exporter->writeStep("initial.vtk"));

      feNorm *norm, *plift, *pdrag, *vlift, *vdrag, *lambdaNorm;
      feCheck(createNorm(norm, VECTOR_L2, {uDomaine}, &sol));
      // feCheck(createNorm(plift, PRESSURE_LIFT_FORCE, {pCylindre}, &sol));
      // feCheck(createNorm(vlift,  VISCOUS_LIFT_FORCE, {uNoSlip  }, &sol));
      // feCheck(createNorm(pdrag, PRESSURE_DRAG_FORCE, {pCylindre}, &sol));
      // feCheck(createNorm(vdrag,  VISCOUS_DRAG_FORCE, {uNoSlip  }, &sol));
      // feCheck(createNorm(lambdaNorm,  L2_ERROR, {lambdaCylindre}, &sol));
      std::vector<feNorm*> norms = {norm};

      feTolerances tol{1e-10, 1e-10, 200};

      TimeIntegrator *solver;

      // if(iAdapt <= nStepsToIncreaseRe) {
        // Start from steady solution at low Re
        // TimeIntegrator *solverSteady;
        // feCheck(createTimeIntegrator(solverSteady, STATIONARY, tol, linearSystem, &numbering, &sol, &mesh, norms, exportData));
        // feInfo("Solving steady state for Re = %f", Re);
        // feCheck(solverSteady->makeStep());
        // for(auto *space : spaces)
        //   space->setDOFInitialization(dofInitialization::PREVIOUS_SOL);
        // uIn->setDOFInitialization(dofInitialization::NODEWISE);
        // container = solverSteady->getSolutionContainer();
        // exit(-1);
      // } else {

        double t0 = theT0;
        double t1 = t0 + dt*nTimeSteps;
        feCheck(createTimeIntegrator(solver, BDF2, tol, linearSystem, &numbering, &sol, &mesh, norms, exportData, t0, t1, nTimeSteps));

        if(iAdapt > 0) {
          solver->setCurrentStep(currentStep);
        }

        // Solve
        for(int ii = 0; ii < nTimeSteps; ++ii) {
          feInfo("Solving unsteady for Re = %f", Re);
          feCheck(solver->makeStep());
        }

        currentStep = solver->getCurrentStep();
        theT0 = sol.getCurrentTime();
        container = solver->getSolutionContainer();
      // }

      // feInfo("Computed pdrag %f (ex)", pdrag->compute());
      // feInfo("Computed plift %f (ey)", plift->compute());
      // feInfo("Computed vdrag %f (ex)", vdrag->computeViscousDrag(D/Re, &recU, &recV));
      // feInfo("Computed total force %f (ex)", lambdaNorm->computeForcesFromLagrangeMultiplier(0));
      // feInfo("Computed total force %f (ey)", lambdaNorm->computeForcesFromLagrangeMultiplier(1));

      // forces_x[iAdapt] = lambdaNorm->computeForcesFromLagrangeMultiplier(0);
      // forces_y[iAdapt] = lambdaNorm->computeForcesFromLagrangeMultiplier(1);
      // ===== Solver ends here ==========================

      if(adapt) {
        // ===== Mesh adaptation starts here ==========================
        // Don't adapt at the last loop, only compute the solution on the last adapted mesh
        if(iAdapt < nAdaptationCycles) {

          bool reconstructAtHighOrderNodes = false;
          feNewRecovery recU(uDomaine, 0, &mesh, &sol, meshFile, "recoveredDerivativesU.msh", reconstructAtHighOrderNodes, false, nullptr, &numbering);
          // feNewRecovery recV(uDomaine, 1, &mesh, &sol, meshFile, "recoveredDerivativesV.msh", reconstructAtHighOrderNodes, false, nullptr, &numbering);
          std::vector<feNewRecovery*> recoveredFields = {&recU};

          feMetricOptions options(meshFile);
          options.nTargetVertices = targetN[iV];
          options.LpNorm = 100.;
          options.insideCallback = inside;
          options.hMin = 1e-10;
          options.hMax = 10.;

          curveTheMesh = false;
            
          if(orderVelocity == 1) {
            options.method = adaptationMethod::ANISO_P1;
            options.enableGradation = true;
            options.gradation = 1.5;
          } else {
            if(curveTheMesh) {
              options.method = adaptationMethod::CURVED_EXTREME_SIZES;
              options.enableGradation = true;
              options.gradation = 1.3;
              // options.eTargetError = 1e-3;
            } else {
              options.method = adaptationMethod::ANISO_PN;
              options.enableGradation = true;
              options.gradation = 1.3;
              options.logSimplexOptions.nThetaPerQuadrant = 20;
              options.logSimplexOptions.tol = 0.02;
              options.logSimplexOptions.uniformErrorCurve = true;
              options.logSimplexOptions.numLoopsUniformErrorCurve = 5;
            }
          }

          bool isBackMeshP2 = false;
          bool setGmshModelToP1 = true; // Set this to false only to test convergence of metric interpolation
          bool curveMMGmesh = true;

          feCheck(mesh.adapt(recoveredFields, options, spaces, essentialSpaces, uDomaine,
            nullptr, nullptr, nullptr, curveTheMesh, isBackMeshP2, setGmshModelToP1, curveMMGmesh, curveToMinimize::DISTORTION));

          // Interpolate solution on adapted mesh
          feMesh2DP1 nextMesh("anisoadapted.msh");
          // feMesh2DP1 nextMesh("../data/cylindre.msh");

          feSpace *uDomaine = nullptr, *pDomaine = nullptr, *uIn = nullptr, *uTop = nullptr, *uBottom = nullptr, *uNoSlip = nullptr, *wDomaine = nullptr;

          if(iAdapt <= nStepsToIncreaseRe)
            feCheck(createFiniteElementSpace(uIn, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Gauche", degreeQuadrature, &inletUniform));
          else
            feCheck(createFiniteElementSpace(uIn, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Gauche", degreeQuadrature, &inletUniform));

          // Cylindre
          // feCheck(createFiniteElementSpace(  uTop, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Haut", degreeQuadrature, &zeroVector));
          // feCheck(createFiniteElementSpace(uBottom, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Bas", degreeQuadrature, &zeroVector));
          // feCheck(createFiniteElementSpace(uNoSlip, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Cylindre", degreeQuadrature, &zeroVector));
          // feCheck(createFiniteElementSpace(uDomaine, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domaine", degreeQuadrature, &inletUniform));
          // feCheck(createFiniteElementSpace(pDomaine, &nextMesh, elementType::LAGRANGE, orderPressure, "P", "Domaine", degreeQuadrature, &zero));
          // // feCheck(createFiniteElementSpace(wDomaine, &nextMesh, elementType::LAGRANGE, orderVelocity, "W", "Domaine", degreeQuadrature, &zero));
          // std::vector<feSpace*> nextSpaces = {uTop, uBottom, uNoSlip, uDomaine, uIn, pDomaine};
          // std::vector<feSpace*> nextEssentialSpaces = {uNoSlip, uIn};
          // uTop->setEssentialComponent(1, true);
          // uBottom->setEssentialComponent(1, true);

          // Poiseuille
          feCheck(createFiniteElementSpace(     uIn, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Gauche", degreeQuadrature, &inletUniform));
          feCheck(createFiniteElementSpace( uDroite, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "Droite", degreeQuadrature, &zeroVector));
          feCheck(createFiniteElementSpace( uNoSlip, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U",  "NoSlip", degreeQuadrature, &zeroVector));
          feCheck(createFiniteElementSpace(uDomaine, &nextMesh, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domaine", degreeQuadrature, &zeroVector));
          feCheck(createFiniteElementSpace(pDomaine, &nextMesh, elementType::LAGRANGE,        orderPressure, "P", "Domaine", degreeQuadrature, &zero));
          std::vector<feSpace*> nextSpaces = {uIn, uDroite, uNoSlip, uDomaine, pDomaine};
          std::vector<feSpace*> nextEssentialSpaces = {uIn, uNoSlip};

          feMetaNumber nextNumbering(&nextMesh, nextSpaces, nextEssentialSpaces);
          feSolution nextSol(nextNumbering.getNbDOFs(), nextSpaces, nextEssentialSpaces);

          // Interpolate previous solution on current adapted mesh
          feCheck(mesh.transfer(&nextMesh, &numbering, &nextNumbering, container, spaces, essentialSpaces, nextSpaces));

          nextSol.setSolFromContainer(container, 0);

          feCheck(exporter->writeStep("afterTransfer.vtk"));

          // Save to final aniso mesh
          std::string cmd = "cp anisoadapted.msh cylinder_finalAniso_" + std::to_string(targetN[iV]) + ".msh";
          system(cmd.data());
        } else {

          if(iV > 0) {

            feInfo("Should not be here unless mesh convergence study");
            exit(-1);
            // Interpolate converged solution from mesh N_i on mesh N_i+1, then compute error

            // Not super efficient, but re-read the previous converged mesh and recreate the structures to project on current mesh
            std::string meshNname = "cylinder_finalAniso_" + std::to_string(targetN[iV-1]) + ".msh";
            feMesh2DP1 meshN(meshNname);

            feSpace *uDomaine = nullptr, *pDomaine = nullptr, *uIn = nullptr, *uTop = nullptr, *uBottom = nullptr, *uNoSlip = nullptr, *wDomaine = nullptr;
            feCheck(createFiniteElementSpace(uIn, &meshN, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Gauche", degreeQuadrature, &inlet));
            feCheck(createFiniteElementSpace(uTop, &meshN, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Haut", degreeQuadrature, &zeroVector));
            feCheck(createFiniteElementSpace(uBottom, &meshN, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Bas", degreeQuadrature, &zeroVector));
            feCheck(createFiniteElementSpace(uNoSlip, &meshN, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Cylindre", degreeQuadrature, &zeroVector));
            feCheck(createFiniteElementSpace(uDomaine, &meshN, elementType::VECTOR_LAGRANGE, orderVelocity, "U", "Domaine", degreeQuadrature, &zeroVector));
            feCheck(createFiniteElementSpace(pDomaine, &meshN, elementType::LAGRANGE, orderPressure, "P", "Domaine", degreeQuadrature, &zero));
            // feCheck(createFiniteElementSpace(wDomaine, &meshN, elementType::LAGRANGE, orderVelocity, "W", "Domaine", degreeQuadrature, &zero));
            std::vector<feSpace*> spacesN = {uIn, uTop, uBottom, uNoSlip, uDomaine, pDomaine};
            std::vector<feSpace*> essentialSpacesN = {uIn, uNoSlip};
            uTop->setEssentialComponent(1, true);
            uBottom->setEssentialComponent(1, true);

            feMetaNumber numberingN(&meshN, spacesN, essentialSpacesN);
            feSolution solN(numberingN.getNbDOFs(), spacesN, essentialSpacesN);
            
            solN.setSolFromContainer(&containerN);

            // Just to debug
            feExporter *exporterN;
            feCheck(createVisualizationExporter(exporterN, VTK, &numberingN, &solN, &meshN, spacesN));
            feCheck(exporterN->writeStep("previousConverged.vtk"));
            
            // Estimate error by projecting previous solution on current mesh (current mesh is mesh2N)
            meshN.transfer(&mesh, &numberingN, &numbering, &containerN, spacesN, essentialSpacesN, spaces);

            // Just to debug
            // Update current solution to plot the interpolated coarse solution on this mesh
            feSolution solN_on_2N;
            solN_on_2N.setSolFromContainer(&containerN, 0);
            feExporter *exporter2N;
            feCheck(createVisualizationExporter(exporter2N, VTK, &numbering, &solN_on_2N, &mesh, spaces));
            feCheck(exporter2N->writeStep("coarseSolutionOnCurrentMesh.vtk"));
            
            // Compute error
            feNorm *norm;
            feCheck(createNorm(norm, L2_ERROR, {uDomaine}, &sol));
            double err = norm->computeLpErrorFromTransferredSolution(2, &solN_on_2N);
            E[iV] = err;
          }

          // Save current final solution to estimate error on finer mesh
          // containerN is now container2N
          containerN = *(solver->getSolutionContainer());
          // Deep copy
          containerN.copy(*(solver->getSolutionContainer()));
        }
        // ===== Mesh adaptation and error estimation end here =================================
      }

      if(Re < Re_final)
        Re += dRe;

      delete uDomaine;
      delete diffU;
      delete divU;
      delete gradP;
      // delete linearSystem;
      // delete exporter;
      // delete solver;
    }

    for(int i = 0; i < nAdaptationCycles; ++i) {
      feInfo("F_x = %+-1.10e - F_y = %+-1.10e", forces_x[i], forces_y[i]);
    }

    if(iV > 0) {
      // Write error to file
      std::ofstream outfile;
      if(iV == 1 && !appendToErrorFile)
        outfile.open("error.txt");
      else 
        outfile.open("error.txt", std::ios_base::app);
        outfile << nV[iV] << "\t"
      << std::setprecision(4) << E[iV] << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"

      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"

      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t" 

      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"
      << std::setprecision(4) << 0. << "\t"

      << std::setprecision(4) << 0. << "\t"

      << std::setprecision(4) << 0. << "\t" << 0. << "\t"
      << std::setprecision(14) <<  0. << "\t"
      << std::setprecision(14) <<  0. << "\t"
      << std::setprecision(6) <<  0. << "\n";
    }
  }

  for(int i = 0; i < nMesh; ++i)
    feInfo("Computed error %+-1.6e", E[i]);

  petscFinalize();
  return 0;
}