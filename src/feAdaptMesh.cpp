#include "feAdaptMesh.h"
#include "feMetricTools.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

feRecovery *activeRecovery;
feSpace *activeIntSpace;
feNumber *activeNumbering;
feSolution *activeSolution;
feFunction *exactSolution;

#if defined(HAVE_GMSH)
double pointwiseErrorCallback(double *x)
{
  std::vector<double> pos(3, 0.);
  pos[0] = x[0];
  pos[1] = x[1];

  // double recovery = activeRecovery->evalDerivative(0, pos);
  double recovery = exactSolution->eval(0, pos);
  // The finite element solution interpolated at x
  double uh = activeIntSpace->interpolateSolution(activeNumbering, activeSolution, pos);

  return recovery - uh;
}

double errorSquaredCallback(double *xa, double *xb, double *xc, double *xab, double *xbc,
                            double *xca)
{
  int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
  std::vector<double> basisFunctions;
  int numComponents, numOrientations;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  std::vector<double> gradBasisFunctions;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "GradLagrange", numComponents,
                                       gradBasisFunctions, numOrientations);

  double F[6] = {f(activeRecovery, xa[0], xa[1]),   f(activeRecovery, xb[0], xb[1]),
                 f(activeRecovery, xc[0], xc[1]),   f(activeRecovery, xab[0], xab[1]),
                 f(activeRecovery, xbc[0], xbc[1]), f(activeRecovery, xca[0], xca[1])};
  double X[6] = {xa[0], xb[0], xc[0], xab[0], xbc[0], xca[0]};
  double Y[6] = {xa[1], xb[1], xc[1], xab[1], xbc[1], xca[1]};
  double e2 = 0;
  for(size_t i = 0; i < weights.size(); i++) {
    double interpolated = 0, x = 0, y = 0, dxdu = 0, dxdv = 0, dydu = 0, dydv = 0;

    for(size_t j = 0; j < 6; j++) {
      x += basisFunctions[6 * i + j] * X[j];
      y += basisFunctions[6 * i + j] * Y[j];
      dxdu += gradBasisFunctions[3 * (6 * i + j) + 0] * X[j];
      dxdv += gradBasisFunctions[3 * (6 * i + j) + 1] * X[j];
      dydu += gradBasisFunctions[3 * (6 * i + j) + 0] * Y[j];
      dydv += gradBasisFunctions[3 * (6 * i + j) + 1] * Y[j];
      // interpolated += basisFunctions[6 * i + j] * F[j];
    }

    double detJ = fabs(dxdu * dydv - dxdv * dydu);

    std::vector<double> pos(3, 0.);
    pos[0] = x;
    pos[1] = y;

    // The recovered solution interpolated at x
    // double recovery = activeRecovery->evalDerivative(0, pos);
    double recovery = exactSolution->eval(0, pos);
    // The finite element solution interpolated at x
    double uh = activeIntSpace->interpolateSolution(activeNumbering, activeSolution, pos);

    e2 += weights[i] * (recovery - uh) * (recovery - uh) * detJ;
  }
  // printf("Coucou\n");
  return e2;
}

void gradErrorSquaredCallback(double *xa, double *xb, double *xc, double *xab, double *xbc,
                              double *xca, std::vector<double> &grad)
{
  int triangleP2 = gmsh::model::mesh::getElementType("Triangle", 2);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triangleP2, "Gauss12", localCoord, weights);
  std::vector<double> basisFunctions;
  int numComponents, numOrientations;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  std::vector<double> gradBasisFunctions;
  gmsh::model::mesh::getBasisFunctions(triangleP2, localCoord, "GradLagrange", numComponents,
                                       gradBasisFunctions, numOrientations);

  double F[6] = {f(activeRecovery, xa[0], xa[1]),   f(activeRecovery, xb[0], xb[1]),
                 f(activeRecovery, xc[0], xc[1]),   f(activeRecovery, xab[0], xab[1]),
                 f(activeRecovery, xbc[0], xbc[1]), f(activeRecovery, xca[0], xca[1])};
  double X[6] = {xa[0], xb[0], xc[0], xab[0], xbc[0], xca[0]};
  double Y[6] = {xa[1], xb[1], xc[1], xab[1], xbc[1], xca[1]};

  grad[0] = 0.0;
  grad[1] = 0.0;

  for(size_t i = 0; i < weights.size(); i++) {
    double interpolated = 0, x = 0, y = 0, dxdu = 0, dxdv = 0, dydu = 0, dydv = 0;

    for(size_t j = 0; j < 6; j++) {
      x += basisFunctions[6 * i + j] * X[j];
      y += basisFunctions[6 * i + j] * Y[j];
      dxdu += gradBasisFunctions[3 * (6 * i + j) + 0] * X[j];
      dxdv += gradBasisFunctions[3 * (6 * i + j) + 1] * X[j];
      dydu += gradBasisFunctions[3 * (6 * i + j) + 0] * Y[j];
      dydv += gradBasisFunctions[3 * (6 * i + j) + 1] * Y[j];
      // interpolated += basisFunctions[6 * i + j] * F[j];
    }

    double detJ = fabs(dxdu * dydv - dxdv * dydu);

    double drdx = dydv / detJ;
    double drdy = -dxdv / detJ;
    double dsdx = -dydu / detJ;
    double dsdy = dxdu / detJ;

    std::vector<double> pos(3, 0.);
    pos[0] = x;
    pos[1] = y;

    double recovery = activeRecovery->evalDerivative(0, pos);
    // The gradient of the recovered solution interpolated at x
    double gradXRecovery = activeRecovery->evalDerivative(1, pos);
    double gradYRecovery = activeRecovery->evalDerivative(2, pos);

    double uh = activeIntSpace->interpolateSolution(activeNumbering, activeSolution, pos);
    // The gradient of the finite element solution interpolated at x
    std::vector<double> gradrs_uh(3, 0.0);
    activeIntSpace->interpolateSolution_gradrs(activeNumbering, activeSolution, pos, gradrs_uh);
    double duhdx = gradrs_uh[0] * drdx + gradrs_uh[1] * dsdx;
    double duhdy = gradrs_uh[0] * drdy + gradrs_uh[1] * dsdy;

    grad[0] += weights[i] * detJ * 2.0 * (recovery - uh) * (gradXRecovery - duhdx);
    grad[1] += weights[i] * detJ * 2.0 * (recovery - uh) * (gradYRecovery - duhdy);
  }
}
#endif

/* Creates an adapted straight-sided anisotropic mesh based on the computed metric field
*/
void createAnisoMesh(feMetric *metric, feMetricOptions metricOptions)
{
  // Write size map
  metric->writeSizeFieldSol2D("sizeMapAniso.sol");
  // Get input and output mesh names
  size_t lastindex = metricOptions.meshName.find_last_of(".");
  std::string root = metricOptions.meshName.substr(0, lastindex);
  std::string meshToAdapt3D = root + ".mesh";
  std::string meshToAdapt2D = root + "2D.mesh";
  // Convert input 3D .mesh to 2D .mesh
  std::string cmd = "mmg3Dto2D " + meshToAdapt3D;
  system(cmd.c_str());
  // Adapt mesh
  cmd = "mmg2d " + meshToAdapt2D + " -sol sizeMapAniso.sol -hgrad 3 -o " + metricOptions.adaptedMeshName;
  system(cmd.c_str());
  // Open mesh
  cmd = "gmsh " + metricOptions.adaptedMeshName + " &";
  system(cmd.c_str());
}

/* Creates a curved mesh based on
    - the geometry stored in the active Gmsh model with name metricOptions.gmshModel
    - the metric field stored as a view in this model
    - the "inside" callback metricOptions.inside, returning true if a point is inside the geometry
*/
void createCurvedMesh(feFunction *solExact, feMetaNumber *metaNumber, feSolution *sol,
                      feSpace *intSpace, feRecovery *recovery, feMetric *metric,
                      feMetricOptions metricOptions)
{
#if defined(HAVE_GMSH)
  // if(metricOptions.isGmshModelReady){
  std::vector<double> pts;
  std::string modelForMetric = metricOptions.gmshModel;
  std::string modelForMesh = "adapt";
  int faceTag = 0;
  activeRecovery = recovery;
  activeIntSpace = intSpace;
  activeNumbering = metaNumber->getNumbering(intSpace->getFieldID());
  activeSolution = sol;
  exactSolution = solExact;

  gmsh::model::add(modelForMesh);
  gmsh::model::setCurrent(modelForMesh);
  system("mmg2d metric.msh -hgrad 3");
  gmsh::merge("metric.o.msh");

  gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);

  // computePointsUsingScaledCrossFieldPlanarP2(
    // modelForMetric.c_str(), modelForMesh.c_str(), metric->getMetricViewTag(), faceTag, pts,
    // errorSquaredCallback, metricOptions.inside, gradErrorSquaredCallback);

  // computePointsUsingScaledCrossFieldPlanarP2(modelForMetric.c_str(), modelForMesh.c_str(),
  // metric->getMetricViewTag(), faceTag, pts, NULL, metricOptions.inside);

  gmsh::write("sol_adapt.msh");
  system("gmsh sol_adapt.msh &");
  // } else{
  // 	printf("In feAdaptMesh : Error - Metric file is not ready to generate a curved mesh.\n");
  // }
#else
  printf("In feAdaptMesh : Error - Gmsh is required to generate curved meshes.\n");
#endif
}