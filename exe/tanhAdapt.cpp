/*

  Curvilinear mesh adaptation.
  Requires Gmsh.

*/

#include "feAPI.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

static double xmin = -1.;
static double xmax =  1.;
static double ymin = -1.;
static double ymax =  1.;
// static double x_0 = -0.2;
// static double y_0 = -0.2;


bool inside(double* x, bool strictly = false){
  if(strictly){
    return (x[0] > xmin) && (x[1] > ymin) && (x[0] < xmax) && (x[1] < ymax);
  } else{
    return (x[0] >= xmin) && (x[1] >= ymin) && (x[0] <= xmax) && (x[1] <= ymax);
  }
}

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];

  // Quarter on [0,1]
  // double a = 10.0;
  // double b = 2.;
  // return 0.5*(1.0 + tanh(a*(x - sin(b*M_PI*y)/4.0)));

  // return atan(10. * (sin(3.*M_PI*y/2.) - 2.*x));

  // r^4
  // double r = sqrt(x*x + y*y);
  // return (x*x + y*y) * (x*x + y*y);

  // x^3
  // return x*x*x + 20*y*y*y;

  // atan (IMR paper)
  return atan(10.*(sin(3.*M_PI*y/2.) - 2.*x));

  // The full atan on [-1,1]
  // double a = 10.0;
  // double b = 1.5;
  // return 0.5*(1.0 + tanh(a*(x/2. - sin(b*M_PI*y)/4.0 )));

  // return cos(2.0*M_PI*x) * sin(2.0*M_PI*y);

  // double b = 5.0;
  // return sqrt((1.+b*b)/(1.+b*b*sin(2.*M_PI*x)*cos(2.*M_PI*y)*sin(2.*M_PI*x)*cos(2.*M_PI*y))) * sin(2.*M_PI*x)*cos(2.*M_PI*y); // flattenedCos
}

double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  // double k = par[0];
  return 0.;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

int main(int argc, char **argv)
{
  static const char *meshFile = "square.msh";
  int order = 2;
  int degreeQuadrature = 10;
  double eTargetError = 1.;
  int nTargetVertices = 100;
  int nLoopsAnisoMesh = 0;
  int nAdaptationCycles = 1;
  double gradation = -1.;

  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&eTargetError, "-e", "--target_error", "Target error on the adapted mesh edges");
  options.addOption(&gradation, "-grad", "--gradation", "Metric gradation (> 1). Gradation is disabled if no gradation is given or if grad < 0.");
  options.addOption(&nTargetVertices, "-n", "--target_vertices", "Target number of vertices in the adapted mesh");
  options.addOption(&nLoopsAnisoMesh, "-nLoop", "--loop_aniso_mesh", "Number of loops to perform to create the aniso mesh");
  options.addOption(&nAdaptationCycles, "-nAdapt", "--nAdaptation", "Number of adaptation loops");
  feCheck(options.parse());

  std::string meshFileTmp = meshFile;

  setVerbose(2);

  double k = 1.0;
  feFunction funSol(fSol, {});
  feFunction funSource(fSource, {k});
  feFunction funZero(fZero, {});

  feSpace *uBord, *uDomaine;
  feRecovery *recovery;

  // Start the adaptation loop
  for(int iAdapt = 0; iAdapt < nAdaptationCycles; ++iAdapt)
  {
    // bool reversed = iAdapt
    feMesh2DP1 mesh(meshFile);
    feCheck(createFiniteElementSpace(uBord, &mesh, 1, LINE, LAGRANGE, order, "U", "Bord", degreeQuadrature, &funSol));
    feCheck(createFiniteElementSpace(uDomaine, &mesh, 2, TRI, LAGRANGE, order, "U", "Domaine", degreeQuadrature, &funSol));
    std::vector<feSpace *> spaces = {uBord, uDomaine};
    std::vector<feSpace *> essentialSpaces = {uBord, uDomaine};

    feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);
    feSolution sol(metaNumber.getNbDOFs(), spaces, essentialSpaces);
    sol.initializeUnknowns(&mesh);
    sol.initializeEssentialBC(&mesh);

    feExporter *exporter;
    feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
    int exportEveryNSteps = 1;
    std::string vtkFileRoot = "tanh";
    std::string vtkFile = vtkFileRoot + "_" + std::to_string(iAdapt+1) + ".vtk";
    feExportData exportData = {exporter, exportEveryNSteps, vtkFileRoot};
    feCheck(exporter->writeStep(vtkFile));

    // Mesh adaptation
    std::string metricMeshName = "recoveredDerivatives.msh";
    std::vector<double> estErreur(2, 0.);

    recovery = new feRecovery(&metaNumber, uDomaine, &mesh, &sol, estErreur, nullptr, meshFile, metricMeshName);

    feMetricOptions metricOptions;
    // if( iAdapt + 1 == nAdaptationCycles){
      // JF's metric for curving
      // metricOptions.computationMethod = 4;
      // metricOptions.eTargetError = eTargetError;
      // metricOptions.nTargetVertices = -1;
    // } else{
    //   // Logsimplex metric for the intermediary aniso meshes
      metricOptions.computationMethod = 3;
      metricOptions.eTargetError = 1.;
      metricOptions.nTargetVertices = nTargetVertices;
    // }
    
    metricOptions.directionFieldFromDerivativesOfOrder = 1;
    metricOptions.polynomialDegree = 2;

    if(gradation < 0){
      metricOptions.enableGradation = false;
    } else{
      metricOptions.enableGradation = true;
      metricOptions.gradation = gradation;
    }

    metricOptions.LpNorm = 2;

    double modelSize = fmax(ymax - ymin, xmax - xmin);
    metricOptions.modelSizeX = xmax;
    metricOptions.modelSizeY = ymax;
    metricOptions.inside = inside;
    metricOptions.hMin = modelSize / 100000;
    metricOptions.hMax = modelSize / 1.;

  #if defined(HAVE_GMSH)  

    // Open current mesh as background mesh
    gmsh::initialize();
    meshFile = meshFileTmp.c_str();
    gmsh::open(meshFile);
    // gmsh::fltk::run();
    gmsh::model::getCurrent(metricOptions.modelForMetric);

    gmsh::write("thegmshModel.msh");
    metricOptions.isGmshModelReady = true;

    metricOptions.modelForMesh = "modelForMesh";

    metricOptions.adaptedMeshName = "adapted_" + std::to_string(iAdapt+1) + ".msh";
    metricOptions.metricMeshNameForMMG = "metricMMG_" + std::to_string(iAdapt+1) + ".msh";
    metricOptions.metricMeshNameForMMG_out = "outMetricMMG_" + std::to_string(iAdapt+1) + ".msh";
    metricOptions.recoveryName = "recovery_" + std::to_string(iAdapt+1) + ".msh";

    feMetric metric(recovery, metricOptions);
    metric.computeMetrics();

    feInfo("Creating curved mesh...");
    createCurvedMesh(nullptr, &metaNumber, &sol, uDomaine, recovery, &metric, metricOptions, 0, nLoopsAnisoMesh, (iAdapt+1) == nAdaptationCycles);
    // createCurvedMesh(nullptr, &metaNumber, &sol, uDomaine, recovery, &metric, metricOptions, 0, nLoopsAnisoMesh, false);

    // Attention : meshFile will be empty once metricOptions is deleted...
    meshFile = metricOptions.adaptedMeshName.c_str();
    meshFileTmp = metricOptions.adaptedMeshName;
    
    gmsh::clear();
    gmsh::finalize();

  #endif

    delete recovery;
    delete uBord;
    delete uDomaine;
  }

  return 0;
}