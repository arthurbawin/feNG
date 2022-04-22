/*

  Curvilinear mesh adaptation.
  Requires Gmsh.

*/

#include "feAPI.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

static double xmin = 0.;
static double xmax = 1.;
static double ymin = 0.;
static double ymax = 1.;
static double x_0 = -0.2;
static double y_0 = -0.2;


bool inside(double* x){
  bool isInsideRectangle = (x[0] >= xmin) && (x[1] >= ymin) && (x[0] <= xmax) && (x[1] <= ymax);
  return isInsideRectangle;
}

double fSol(const double t, const std::vector<double> &pos, const std::vector<double> &par) {
  double x = pos[0];
  double y = pos[1];
  // return atan(20.0*x - 10.0*sin((3.0*M_PI*y)/2.0));
  // return exp(-10.*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)));
  double rSquared = (x-x_0)*(x-x_0) + (y-y_0)*(y-y_0);
  return pow(rSquared, 2);
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
  feRecovery *recU = new feRecovery(&metaNumber, uDomaine, &mesh, &sol, estErreur, funZero, meshFile, metricMeshName);

  feMetricOptions metricOptions;
  metricOptions.computationMethod = 1;
  metricOptions.polynomialDegree = 2;
  metricOptions.eTargetError = eTargetError;
  metricOptions.nTargetVertices = nTargetVertices;
  metricOptions.nPhi = 501;

  double modelSize = fmax(ymax - ymin, xmax - xmin);
  metricOptions.modelSizeX = xmax;
  metricOptions.modelSizeY = ymax;
  metricOptions.inside = inside;

  metricOptions.hMin = modelSize / 100000;
  metricOptions.hMax = modelSize;

#if defined(HAVE_GMSH)
  // // Create a gmsh model with the geometry
  // metricOptions.gmshModel = "myGeometry";
  // gmsh::initialize();
  // gmsh::model::add(metricOptions.gmshModel);
  // double xmax = metricOptions.modelSizeX;
  // double ymax = metricOptions.modelSizeY;
  // double lc = modelSize/50.0;
  // gmsh::model::occ::addPoint(0, 0, 0, lc, 1);
  // gmsh::model::occ::addPoint(xmax, 0, 0, lc, 2);
  // gmsh::model::occ::addPoint(xmax, ymax, 0, lc, 3);
  // gmsh::model::occ::addPoint(0, ymax, 0, lc, 4);
  // int l1 = gmsh::model::occ::addLine(1, 2);
  // int l2 = gmsh::model::occ::addLine(2, 3);
  // int l3 = gmsh::model::occ::addLine(3, 4);
  // int l4 = gmsh::model::occ::addLine(4, 1);
  // int boundary = gmsh::model::occ::addCurveLoop({l1, l2, l3, l4});
  // gmsh::model::occ::addPlaneSurface({boundary}, 1);
  // gmsh::model::occ::synchronize();
  // gmsh::option::setNumber("Mesh.MeshSizeMin", modelSize/50);
  // gmsh::option::setNumber("Mesh.MeshSizeMax", modelSize/50);
  // gmsh::model::mesh::generate();
  // // gmsh::fltk::run();
  // metricOptions.isGmshModelReady = true;

  // metricOptions.directionFieldFromDerivativesOfOrder = 1;

  feMetric metric(recU, metricOptions);
  metric.computeMetrics();

  metricOptions.meshName = meshFile;
  metricOptions.adaptedMeshName = "foo.mesh";

  feInfo("Creating aniso mesh...");
  createAnisoMesh(&metric, metricOptions);
  // createCurvedMesh(funSol, &metaNumber, &sol, uDomaine, recU, &metric, metricOptions);
#endif

  delete exporter;
  delete uBord;
  delete uDomaine;
  delete funSol;
  delete funSource;

  return 0;
}