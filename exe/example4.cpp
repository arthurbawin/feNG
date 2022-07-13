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
static double x_0 = -0.2;
static double y_0 = -0.2;


bool inside(double* x){
  bool isInsideRectangle = (x[0] >= xmin) && (x[1] >= ymin) && (x[0] <= xmax) && (x[1] <= ymax);
  return isInsideRectangle;
}

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];

  // Quarter on [0,1]
  // double a = 10.0;
  // double b = 2.;
  // return 0.5*(1.0 + tanh(a*(x - sin(b*M_PI*y)/4.0)));

  // The full atan on [-1,1]
  double a = 10.0;
  double b = 1.5;
  return 0.5*(1.0 + tanh(a*(x/2. - sin(b*M_PI*y)/4.0)));

  // return cos(2.0*M_PI*x) * sin(2.0*M_PI*y);

  // double b = 5.0;
  // return sqrt((1.+b*b)/(1.+b*b*sin(2.*M_PI*x)*cos(2.*M_PI*y)*sin(2.*M_PI*x)*cos(2.*M_PI*y))) * sin(2.*M_PI*x)*cos(2.*M_PI*y); // flattenedCos
}

double fSource(const double t, const std::vector<double> &pos, const std::vector<double> &par)
{
  double k = par[0];
  return 0.;
}

double fZero(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  return 0.0;
}

int main(int argc, char **argv)
{
  const char *meshFile = "square.msh";
  int order = 2;
  int degreeQuadrature = 10;
  double eTargetError = 1e-2;
  int nTargetVertices = 1000;

  feOptionsParser options(argc, argv);
  options.addOption(&meshFile, "-m", "--mesh", "Mesh file");
  options.addOption(&eTargetError, "-e", "--target_error", "Target error on the adapted mesh edges");
  options.addOption(&nTargetVertices, "-n", "--target_vertices", "Target number of vertices in the adapted mesh");
  feCheck(options.parse());

  setVerbose(2);

  feMesh2DP1 mesh(meshFile);

  double k = 1.0;
  feFunction *funSol = new feFunction(fSol, {});
  feFunction *funSource = new feFunction(fSource, {k});
  feFunction *funZero = new feFunction(fZero, {});

  feSpace *uBord, *uDomaine;
  feCheck(createFiniteElementSpace(uBord, &mesh, 1, LINE, order, "U", "Bord", degreeQuadrature, funSol));
  feCheck(createFiniteElementSpace(uDomaine, &mesh, 2, TRI, order, "U", "Domaine", degreeQuadrature, funSol));

  std::vector<feSpace *> spaces = {uBord, uDomaine};
  std::vector<feSpace *> essentialSpaces = {uBord, uDomaine};

  feMetaNumber metaNumber(&mesh, spaces, essentialSpaces);

  feSolution sol(&mesh, spaces, essentialSpaces, &metaNumber);
  sol.initializeUnknowns(&mesh, &metaNumber);
  sol.initializeEssentialBC(&mesh, &metaNumber);

  feExporter *exporter;
  feCheck(createVisualizationExporter(exporter, VTK, &metaNumber, &sol, &mesh, spaces));
  std::string vtkFileRoot = "root";
  feCheck(exporter->writeStep(vtkFileRoot + "Analytic.vtk"));

  // Mesh adaptation
  std::string metricMeshName = "recoveredDerivatives.msh";
  std::vector<double> estErreur(2, 0.);

  // feRecovery *recU = new feRecovery(&metaNumber, uDomaine, &mesh, &sol, estErreur, funZero, meshFile, metricMeshName);
  // recU->writeRecovery("recoveryFileFullAtan2.txt");
  feRecovery *recU = new feRecovery(uDomaine, &mesh, "recoveryFileFullAtan.txt");
  // recU->writeRecovery("recoveryFile_compare.txt");

  feMetricOptions metricOptions;
  metricOptions.computationMethod = 2;
  metricOptions.polynomialDegree = 2;
  metricOptions.eTargetError = eTargetError;
  metricOptions.nTargetVertices = nTargetVertices;
  metricOptions.nPhi = 501;
  metricOptions.LpNorm = 2;

  double modelSize = fmax(ymax - ymin, xmax - xmin);
  metricOptions.modelSizeX = xmax;
  metricOptions.modelSizeY = ymax;
  metricOptions.inside = inside;

  metricOptions.hMin = modelSize / 100000;
  metricOptions.hMax = modelSize / 1.;

  metricOptions.directionFieldFromDerivativesOfOrder = 1;

#if defined(HAVE_GMSH)  

  for(int iCycle = 0; iCycle < 2; ++iCycle){

    if(iCycle == 0){
      // Create a mostly structured gmsh model with the geometry for the first adaptation
      metricOptions.modelForMetric = "myGeometry";
      gmsh::initialize();
      gmsh::model::add(metricOptions.modelForMetric);
      double xmax = metricOptions.modelSizeX;
      double ymax = metricOptions.modelSizeY;
      double lc = modelSize/40.0;
      gmsh::model::occ::addPoint(xmin, ymin, 0, lc, 1);
      gmsh::model::occ::addPoint(xmax, ymin, 0, lc, 2);
      gmsh::model::occ::addPoint(xmax, ymax, 0, lc, 3);
      gmsh::model::occ::addPoint(xmin, ymax, 0, lc, 4);
      int l1 = gmsh::model::occ::addLine(1, 2);
      int l2 = gmsh::model::occ::addLine(2, 3);
      int l3 = gmsh::model::occ::addLine(3, 4);
      int l4 = gmsh::model::occ::addLine(4, 1);
      int boundary = gmsh::model::occ::addCurveLoop({l1, l2, l3, l4});
      gmsh::model::occ::addPlaneSurface({boundary}, 1);
      gmsh::model::occ::synchronize();
      gmsh::option::setNumber("Mesh.MeshSizeMin", modelSize/40.0);
      gmsh::option::setNumber("Mesh.MeshSizeMax", modelSize/40.0);
      gmsh::model::mesh::generate();
      gmsh::write("thegmshModel.msh");
      // gmsh::fltk::run();
      // gmsh::merge(meshFile);
      metricOptions.isGmshModelReady = true;
    } else{
      // Merge the previous curved mesh and compute metric field on this mesh
      // gmsh::model::add(metricOptions.modelForMetric);
      gmsh::initialize();
      gmsh::open("adapted_" + std::to_string(iCycle) + ".msh");
      gmsh::model::getCurrent(metricOptions.modelForMetric);
      // feInfo("New cycle - Showing the merged previous mesh (current model name is %s)", metricOptions.modelForMetric.c_str());
      gmsh::model::mesh::setOrder(1);
      gmsh::fltk::run();
      gmsh::write("thegmshModel.msh");
      metricOptions.isGmshModelReady = true;
    }

    metricOptions.modelForMesh = "modelForMesh";

    // metricOptions.meshName = meshFile;
    metricOptions.adaptedMeshName = "adapted_" + std::to_string(iCycle+1) + ".msh";
    metricOptions.metricMeshNameForMMG = "metricMMG_" + std::to_string(iCycle+1) + ".msh";
    metricOptions.metricMeshNameForMMG_out = "outMetricMMG_" + std::to_string(iCycle+1) + ".msh";
    metricOptions.recoveryName = "recovery_" + std::to_string(iCycle+1) + ".msh";

    feMetric metric(recU, metricOptions);
    metric.computeMetrics();

    feInfo("Creating curved mesh...");
    // createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions, iCycle != 3);
    // createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions, iCycle != 2);
    // createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions, iCycle != 1);

    if(iCycle == 0) 
      createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions, 0);
    else
      createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions, 10);

    gmsh::clear();
    gmsh::finalize();
  }

  
#endif

  delete exporter;
  delete uBord;
  delete uDomaine;
  delete funSol;
  delete funSource;

  return 0;
}