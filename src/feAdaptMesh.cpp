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
  double uh = activeIntSpace->interpolateField(activeSolution, pos);

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

  double F[6] = {f(activeRecovery, xa),   f(activeRecovery, xb),
                 f(activeRecovery, xc),   f(activeRecovery, xab),
                 f(activeRecovery, xbc), f(activeRecovery, xca)};
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
    double uh = activeIntSpace->interpolateField(activeSolution, pos);

    e2 += weights[i] * (recovery - uh) * (recovery - uh) * detJ;
  }
  // printf("Coucou\n");
  return e2;
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
  cmd = "mmg2d " + meshToAdapt2D + " -sol sizeMapAniso.sol -hgrad 3 -o " +
        metricOptions.adaptedMeshName;
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
                      feMetricOptions &metricOptions, int onlyGenerateVertices, 
                      int nLoopsAnisoMesh, bool curve)
{
#if defined(HAVE_GMSH)
  // if(metricOptions.isGmshModelReady){
  std::vector<double> pts;
  int faceTag = 0;
  activeRecovery = recovery;
  activeIntSpace = intSpace;
  activeNumbering = metaNumber->getNumbering(intSpace->getFieldID());
  activeSolution = sol;
  exactSolution = solExact;

  gmsh::model::add(metricOptions.modelForMesh);
  gmsh::model::setCurrent(metricOptions.modelForMesh);

  // Aniso mesh with MMG (used to get the boundary vertices only)
  // std::string cmd = "mmg2d " + metricOptions.metricMeshNameForMMG + " -hgrad 3 -o " + metricOptions.metricMeshNameForMMG_out;
  // std::string cmd1 = "mmg2d " + metricOptions.metricMeshNameForMMG + " -hgrad 10 -o tmp.mesh";
  std::string cmd1 = "mmg2d " + metricOptions.metricMeshNameForMMG + " -hgrad -1 -o tmp.mesh";
  system(cmd1.c_str());
  std::string cmd2 = "gmsh tmp.mesh -o " + metricOptions.metricMeshNameForMMG_out + " -0";
  system(cmd2.c_str());
  std::string cmd3 = "gmsh " + metricOptions.metricMeshNameForMMG_out + " &";
  // system(cmd3.c_str());

  // Loop a few times
  for(int i = 0; i < nLoopsAnisoMesh; ++i){
    gmsh::open(metricOptions.metricMeshNameForMMG_out);
    gmsh::model::getCurrent(metricOptions.modelForMetric);
    metric->setGmshMetricModel(metricOptions.modelForMetric);
    // feInfo("Showing gmsh models at iter %d", i);
    // gmsh::fltk::run();
    metric->computeMetrics();
    system(cmd1.c_str());
    system(cmd2.c_str());
    // system(cmd3.c_str());
  };

  // system(cmd3.c_str());

  gmsh::clear();
  gmsh::open(metricOptions.metricMeshNameForMMG_out);
  // gmsh::open("thegmshModel.msh");
  gmsh::model::getCurrent(metricOptions.modelForMesh);

  gmsh::write("beforeVersion.msh");
  gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
  gmsh::write("afterVersion.msh");

  // Determine if the mesh is reversed (MMG seems to reverse the mesh sometimes, unless I'm doing something wrong)
  // Get quadrature rule and interpolation functions on the adapted mesh
  int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  std::vector<double> localCoord;
  std::vector<double> weights;
  gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss4", localCoord, weights);
  // Get the jacobians
  std::vector<double> jac, det, points;
  gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);

  // for(int i = 0; i < 20; ++i)
  //   feInfo("det1 = %f", det[i]);

  // Check the first determinant
  if(det[0] < 0){
    feInfo("Adapted aniso mesh is reversed : reversing the numbering");
    gmsh::model::mesh::reverse();
  }

  gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);
  // for(int i = 0; i < 20; ++i)
  //   feInfo("det2 = %f", det[i]);

  // Get geometric entities and assign physical entities (should be improved)
  gmsh::vectorpair dimTags;
  gmsh::model::getEntities(dimTags);
  std::vector<int> entities1D, entities2D;
  for(auto p : dimTags){
    std::cout<<p.first<<" - "<<p.second<<std::endl;
    if(p.first == 1)
      entities1D.push_back(p.second);
    if(p.first == 2)
      entities2D.push_back(p.second);
  }

  gmsh::write("beforePhysical.msh");
  // Add physical entities
  gmsh::model::addPhysicalGroup(1, entities1D, 1);
  gmsh::model::setPhysicalName(1, 1, "Bord");
  gmsh::model::addPhysicalGroup(2, entities2D, 2);
  gmsh::model::setPhysicalName(2, 2, "Domaine");
  gmsh::write("afterPhysical.msh");

  // For some reason (?) the physical tags added above can be negative, although
  // the mesh is numbered in counterclockwise orientation.
  // If it's the case, reverse the mesh and keep the negative physical tags.
  std::vector<int> physicalTags;
  bool atLeastOnePositive = false, atLeastOneNegative = false;
  for(auto e : entities1D){
    gmsh::model::getPhysicalGroupsForEntity(1, e, physicalTags);
    if(physicalTags[0] >= 0){
      atLeastOnePositive = true;
    } else{
      atLeastOneNegative = true;
    }
  }
  for(auto e : entities2D){
    gmsh::model::getPhysicalGroupsForEntity(2, e, physicalTags);
    if(physicalTags[0] >= 0){
      atLeastOnePositive = true;
    } else{
      atLeastOneNegative = true;
    }
  }
  if(atLeastOnePositive && atLeastOneNegative){
    // If this happens then I'm very confused
    feWarning("Some physical tags added to the MMG mesh are positive, whereas some are negative...");
    exit(-1);
  }

  if(atLeastOneNegative){
    // All are negative : reverse the mesh
    feInfo("Physical tags added to the MMG mesh are negative but numbering is positive :"
      " reversing the numbering to match the signs of the physical tags");
    gmsh::model::mesh::reverse();
  }
  gmsh::write("afterPhysicalCheck.msh");
  
  gmsh::write(metricOptions.adaptedMeshName);

  // Curve after a few aniso adaptations
  if(curve){
    gmsh::open(metricOptions.metricMeshNameForMMG);
    gmsh::model::getCurrent(metricOptions.modelForMetric);

    gmsh::model::setCurrent(metricOptions.modelForMetric);
    gmsh::model::setCurrent(metricOptions.modelForMesh);

    gmsh::model::getEntities(dimTags, 2);
    if(dimTags.size() > 1){
      feWarning("Gmsh model has more than one surface");
    } else{
      faceTag = dimTags[0].second;
    }

    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    feInfo("There are %d views in the gmsh model : ", viewTags.size());
    for(auto val : viewTags){
      feInfo("View %d with index %d", val, gmsh::view::getIndex(val));
    }

    metric->setMetricViewTag(viewTags[0]);

    computePointsUsingScaledCrossFieldPlanarP2(
      metricOptions.modelForMetric.c_str(),
      metricOptions.modelForMesh.c_str(),
      metric->getMetricViewTag(), 
      faceTag,
      pts,
      errorSquaredCallback,
      metricOptions.inside,
      nullptr,
      onlyGenerateVertices,
      evaluateFieldFromRecoveryCallback,
      (void *) recovery,
      interpolateMetricP1WithDerivativesWrapper,
      interpolateMetricP1Wrapper,
      interpolateMetricAndDerivativeOnP2EdgeWrapper,
      interpolateMetricP1Wrapper1D,
      interpolateMetricAndDerivativeOnP2EdgeWrapper1D,
      (void *) metric);

    gmsh::write(metricOptions.adaptedMeshName);
  }

#else
  printf("In feAdaptMesh : Error - Gmsh is required to generate curved meshes.\n");
#endif
}