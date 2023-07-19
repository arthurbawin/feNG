#include "feMesh.h"
#include "feNorm.h"
#include "feMetric.h"
#include "feMetricTools.h"
#include "feMinimizeInterpolationError.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#include "curvedMesh_structs.h"
#endif

extern int FE_VERBOSE;

int gmshWasInitialized = false;

inline feStatus checkMMGcall(std::string &command)
{
  // From mmg/src/mmg2d/libmmgtypes.h
  int MMG_SUCCESS = 0;

  auto status = system(command.data());
  if(status < 0) {
    return feErrorMsg(FE_STATUS_ERROR, "System call to MMG2D failed with system status %d", status);
  } else {
    if(WIFEXITED(status)) {
      if(WEXITSTATUS(status) == MMG_SUCCESS)
        // Process ran normally and MMG exited normally
        return FE_STATUS_OK;
      else
        // Process ran normally, but MMG failed to adapt for some reason
        return feErrorMsg(FE_STATUS_ERROR, "Adaptation with MMG2D failed with exit code %d",
                          WIFEXITED(status));
    } else
      return feErrorMsg(FE_STATUS_ERROR, "System call to MMG2D failed");
  }
}

feStatus feMesh2DP1::adapt(std::vector<feNewRecovery*> recoveredFields, feMetricOptions &options,
                           const std::vector<feSpace *> &spaces,
                           const std::vector<feSpace *> &essentialSpaces,
                           feSpace *spaceForAdaptation,
                           feSolution *discreteSolution,
                           feFunction *exactSolution,
                           feVectorFunction *exactGradient,
                           bool curve, bool isBackmeshP2, bool setGmshModelToP1,
                           bool curveMMGmesh, curveToMinimize target, feRecovery *oldRecovery)
{

#if defined(HAVE_GMSH)

  activeMesh = this;
  activeRecovery = recoveredFields[0];

  // Step 1: Set the gmsh background mesh on which the metric tensors are computed
  if(!gmshWasInitialized) {

    // Get max threads *before* initializing gmsh
    #if defined(HAVE_OMP)
    int maxNumThreads = omp_get_max_threads();
    #endif

    gmsh::initialize();

    #if defined(HAVE_OMP)
    gmsh::option::setNumber("General.NumThreads", maxNumThreads);
    #endif

    gmshWasInitialized = true;
  }

  if(FE_VERBOSE == VERBOSE_NONE) gmsh::option::setNumber("General.Verbosity", 2);

  gmsh::open(options.backgroundMeshfile);

  int elementOrder = getGeometricInterpolantDegree(recoveredFields[0]->_cnc->getInterpolant());
  if(setGmshModelToP1 && !curve && elementOrder > 1)
  {
    // Set the mesh to P1.
    // This is because MMG only takes a P1 mesh for aniso adaptation.
    // The map from gmsh P1 nodeTags to the feMesh P2 vertices is 
    // built when computing the metrics in feMetric::createVertex2NodeMap
    // If we curve, the background mesh should be P2, so dont change order.

    // For convergence tests on metric interpolation we still want a P2
    // mesh though
    gmsh::model::mesh::setOrder(1);
  }

  gmsh::model::getCurrent(options.modelForMetric);
  options.isGmshModelReady = true;

  feMetric metricField(recoveredFields, options);

  // Step 2: Create aniso mesh
  // Save directly to .msh to preserve Physical Entities
  std::string cmd1 = "mmg2d " + options.mmgInputMeshfile + " -hgrad -1 -o " +
  options.mmgOutputMeshfile;
  if(FE_VERBOSE == VERBOSE_NONE) {
    // Write MMG console outputs to logMMG.txt
    cmd1 += " > logMMG.txt";
  }

  // Assign old feRecovery structure (for curved adaptation, temporary)
  metricField.setRecovery(oldRecovery);

  // The back feMesh can be P2, but we need to write a P1 view in the gmsh model
  // to give mmg for aniso adaptation. So what decides on the number of vertices
  // when metrics are computed is whether we curve (and dont have to give a mesh to mmg)
  // and not whether isP2backmesh is true
  metricField._nVerticesPerElmOnBackmesh = curve ? 6 : 3;

  if(curve && !isBackmeshP2) {
    return feErrorMsg(FE_STATUS_ERROR, "Backmesh must be P2 when curving the mesh.");
  }

  // Compute the metric tensor field and write metrics to options.mmgInputMeshfile
  tic();
  feStatus s = metricField.computeMetrics();
  options.userValue = metricField._options.userValue;

  // To test interpolation error only
  // return FE_STATUS_OK;

  feInfoCond(FE_VERBOSE > 0, "\t\tComputed metric tensors in %f s", toc());
  if(s != FE_STATUS_OK) {
    gmsh::finalize();
    return s;
  }

  if(!curve){
    // Create aniso mesh
    s = checkMMGcall(cmd1);
    if(s != FE_STATUS_OK) {
      gmsh::finalize();
      return s;
    }

    // Open adapted mesh (DEPRECATED: compute the next metric field on this mesh)
    gmsh::open(options.mmgOutputMeshfile);
    gmsh::model::getCurrent(options.modelForMetric);
    metricField.setGmshMetricModel(options.modelForMetric);

    // Step 3: Check physical entities
    gmsh::vectorpair physicalGroups;
    gmsh::model::getPhysicalGroups(physicalGroups);

    std::map<std::pair<int, int>, std::vector<int> > entitiesForPhysical;

    // Get the entities for existing physical groups
    for(auto pair : physicalGroups) {
      std::vector<int> entities;
      gmsh::model::getEntitiesForPhysicalGroup(pair.first, pair.second, entities);
      entitiesForPhysical[pair] = entities;
    }

    // Remove all the physical groups
    gmsh::model::removePhysicalGroups();

    // Re-add stored physical groups
    for(auto pair : _physicalEntitiesDescription) {
      int dim = pair.first.first;
      int tag = pair.first.second;
      std::string name = pair.second;

      bool OK = false;
      for(auto p : physicalGroups) {
        if(p.first == dim && p.second == tag) {
          gmsh::model::addPhysicalGroup(dim, entitiesForPhysical[p], tag);
          gmsh::model::setPhysicalName(dim, tag, name);
          OK = true;
        }
      }
      if(!OK) {
        return feErrorMsg(FE_STATUS_ERROR,
                          "Physical Entity \"%s\" with (dim,tag) = (%d,%d)"
                          " could not be reassigned after mesh adaptation :/",
                          name.data(), dim, tag);
      }
    }

    // gmsh::model::mesh::reverse();

    if(setGmshModelToP1 && !curve && elementOrder > 1)
    {
      // Reset the mesh to Pn
      gmsh::model::mesh::setOrder(elementOrder);
    }

    // Write adapted anisotropic straight mesh with correct Physical Entities
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::write("aniso" + options.adaptedMeshName);
  }

  // Step 4: Interpolate solution on new mesh (TODO)
  // feMesh2DP1 newMesh("afterPhysical.msh");
  // feMetaNumber numbering(&newMesh, spaces, essentialSpaces);

  // Curve after a few aniso adaptations
  if(curve) {

    // Tools to compute interpolation error while curving
    feNorm *norm;
    createNorm(norm, L2_ERROR, {spaceForAdaptation}, discreteSolution, exactSolution, exactGradient);
    norm->setRecovery(recoveredFields[0]);
    activeNorm = norm;
    activeConnectivity = activeMesh->getCncGeoByName("Domaine");
    activeIntSpace = spaceForAdaptation;
    activeSolution = discreteSolution;
    activeExactSolution = exactSolution;
    // activeExactSolution = nullptr;
    activeExactSolutionGradient = exactGradient;
    // activeExactSolutionGradient = nullptr;

    if(target == curveToMinimize::INTERPOLATION_ERROR)
      computeInterpolationErrorOnEachElement();

    if(elementOrder > 1)
    {
      // GFace2PolyMesh only takes P1 meshes
      gmsh::model::mesh::setOrder(1);
    }

    int faceTag;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, 2);
    if(dimTags.size() > 1) {
      return feErrorMsg(FE_STATUS_ERROR, "Gmsh model has more than one surface");
    } else {
      faceTag = dimTags[0].second;
    }

    std::vector<int> viewTags;
    gmsh::view::getTags(viewTags);
    metricField.setMetricViewTag(viewTags[0]);

    // Check there is an "inside" callback
    if(options.insideCallback == nullptr){
      return feErrorMsg(FE_STATUS_ERROR,
        "Provide an \"inside\" callback to peform curved mesh adaptation");
    }

    curvedMeshOptions meshOptions;

    // modelForMetric, modelForMesh and VIEW_TAG are deprecated
    meshOptions.modelForMetric = options.modelForMetric.data();
    meshOptions.modelForMesh = options.modelForMetric.data();
    meshOptions.VIEW_TAG = metricField.getMetricViewTag();

    meshOptions.faceTag = faceTag;
    meshOptions.inside = options.insideCallback;
    meshOptions.curveMMGmesh = curveMMGmesh;
    meshOptions.target = target;
    meshOptions.computeInterpolationError = computeInterpolationError;
    meshOptions.computeInterpolationErrorGradient = computeInterpolationErrorGradient;
    meshOptions.computeInterpolationErrorCallback_EdgesAndVertices = computeInterpolationErrorCallback_EdgesAndVertices;
    meshOptions.applyCurvatureToFeMesh = applyCurvatureToFeMesh;
    meshOptions.getMidnodeTags = getMidnodeTags;
    meshOptions.getPolyMeshVertexTags = getPolyMeshVertexTags;
    meshOptions.metricUserPointer = (void *) &metricField;
    meshOptions.interpolateMetricP1Callback = interpolateMetricP1Callback;
    meshOptions.interpolateMetricP2Callback = interpolateMetricP2Callback;
    meshOptions.interpolateMetricP2CallbackLog = interpolateMetricP2CallbackLog;

    tic();
    createCurvedMesh(meshOptions);
    feInfoCond(FE_VERBOSE > 0, "\t\tGenerated curved mesh in %f s", toc());

    // Metric interpolation error
    // Or
    // Final error after curving when minimizing interpolation error
    options.userValue = meshOptions.userValue;

    gmsh::model::setCurrent(options.modelForMetric.data());
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::write("curvedMesh.msh");
  }

#endif
  return FE_STATUS_OK;
}

/* Creates a curved mesh based on
    - the geometry stored in the active Gmsh model with name metricOptions.gmshModel
    - the metric field stored as a view in this model
    - the "inside" callback metricOptions.insideCallback, returning true if a point is inside the geometry
*/
// void createCurvedMesh(feFunction *solExact, feMetaNumber *metaNumber, feSolution *sol,
//                       feSpace *intSpace, feRecovery *recovery, feMetric *metric,
//                       feMetricOptions &metricOptions, int onlyGenerateVertices, int nLoopsAnisoMesh,
//                       bool curve)
// {
// #if defined(HAVE_GMSH)
//   // if(metricOptions.isGmshModelReady){
//   std::vector<double> pts;
//   int faceTag = 0;
//   activeRecovery = recovery;
//   activeIntSpace = intSpace;
//   activeNumbering = metaNumber->getNumbering(intSpace->getFieldID());
//   activeSolution = sol;
//   exactSolution = solExact;

//   gmsh::model::add(metricOptions.modelForMesh);
//   gmsh::model::setCurrent(metricOptions.modelForMesh);

//   // Aniso mesh with MMG (used to get the boundary vertices only)
//   // std::string cmd = "mmg2d " + metricOptions.mmgInputMeshfile + " -hgrad 3 -o " +
//   // metricOptions.mmgOutputMeshfile; std::string cmd1 = "mmg2d " +
//   // metricOptions.mmgInputMeshfile + " -hgrad 10 -o tmp.mesh";
//   std::string cmd1 = "mmg2d " + metricOptions.mmgInputMeshfile + " -hgrad -1 -o tmp.mesh";
//   system(cmd1.c_str());
//   std::string cmd2 = "gmsh tmp.mesh -o " + metricOptions.mmgOutputMeshfile + " -0";
//   system(cmd2.c_str());
//   std::string cmd3 = "gmsh " + metricOptions.mmgOutputMeshfile + " &";
//   // system(cmd3.c_str());

//   // Loop a few times
//   for(int i = 0; i < nLoopsAnisoMesh; ++i) {
//     gmsh::open(metricOptions.mmgOutputMeshfile);
//     gmsh::model::getCurrent(metricOptions.modelForMetric);
//     metric->setGmshMetricModel(metricOptions.modelForMetric);
//     // feInfo("Showing gmsh models at iter %d", i);
//     // gmsh::fltk::run();
//     metric->computeMetrics();
//     system(cmd1.c_str());
//     system(cmd2.c_str());
//     // system(cmd3.c_str());
//   };

//   // system(cmd3.c_str());

//   gmsh::clear();
//   gmsh::open(metricOptions.mmgOutputMeshfile);
//   // gmsh::open("thegmshModel.msh");
//   gmsh::model::getCurrent(metricOptions.modelForMesh);

//   gmsh::write("beforeVersion.msh");
//   gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
//   gmsh::write("afterVersion.msh");

//   // Determine if the mesh is reversed (MMG seems to reverse the mesh sometimes, unless I'm doing
//   // something wrong) Get quadrature rule and interpolation functions on the adapted mesh
//   int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
//   std::vector<double> localCoord;
//   std::vector<double> weights;
//   gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss4", localCoord, weights);
//   // Get the jacobians
//   std::vector<double> jac, det, points;
//   gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);

//   // for(int i = 0; i < 20; ++i)
//   //   feInfo("det1 = %f", det[i]);

//   // Check the first determinant
//   if(det[0] < 0) {
//     feInfo("Adapted aniso mesh is reversed : reversing the numbering");
//     gmsh::model::mesh::reverse();
//   }

//   gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, points);
//   // for(int i = 0; i < 20; ++i)
//   //   feInfo("det2 = %f", det[i]);

//   // Get geometric entities and assign physical entities (should be improved)
//   gmsh::vectorpair dimTags;
//   gmsh::model::getEntities(dimTags);
//   std::vector<int> entities1D, entities2D;
//   for(auto p : dimTags) {
//     std::cout << p.first << " - " << p.second << std::endl;
//     if(p.first == 1) entities1D.push_back(p.second);
//     if(p.first == 2) entities2D.push_back(p.second);
//   }

//   gmsh::write("beforePhysical.msh");
//   // Add physical entities
//   gmsh::model::addPhysicalGroup(1, entities1D, 1);
//   gmsh::model::setPhysicalName(1, 1, "Bord");
//   gmsh::model::addPhysicalGroup(2, entities2D, 2);
//   gmsh::model::setPhysicalName(2, 2, "Domaine");
//   gmsh::write("afterPhysical.msh");

//   // For some reason (?) the physical tags added above can be negative, although
//   // the mesh is numbered in counterclockwise orientation.
//   // If it's the case, reverse the mesh and keep the negative physical tags.
//   std::vector<int> physicalTags;
//   bool atLeastOnePositive = false, atLeastOneNegative = false;
//   for(auto e : entities1D) {
//     gmsh::model::getPhysicalGroupsForEntity(1, e, physicalTags);
//     if(physicalTags[0] >= 0) {
//       atLeastOnePositive = true;
//     } else {
//       atLeastOneNegative = true;
//     }
//   }
//   for(auto e : entities2D) {
//     gmsh::model::getPhysicalGroupsForEntity(2, e, physicalTags);
//     if(physicalTags[0] >= 0) {
//       atLeastOnePositive = true;
//     } else {
//       atLeastOneNegative = true;
//     }
//   }
//   if(atLeastOnePositive && atLeastOneNegative) {
//     // If this happens then I'm very confused
//     feWarning(
//       "Some physical tags added to the MMG mesh are positive, whereas some are negative...");
//     exit(-1);
//   }

//   if(atLeastOneNegative) {
//     // All are negative : reverse the mesh
//     feInfo("Physical tags added to the MMG mesh are negative but numbering is positive :"
//            " reversing the numbering to match the signs of the physical tags");
//     gmsh::model::mesh::reverse();
//     exit(-1);
//   }
//   gmsh::write("afterPhysicalCheck.msh");

//   gmsh::write(metricOptions.adaptedMeshName);

//   // Curve after a few aniso adaptations
//   if(curve) {
//     gmsh::open(metricOptions.mmgInputMeshfile);
//     gmsh::model::getCurrent(metricOptions.modelForMetric);

//     gmsh::model::setCurrent(metricOptions.modelForMetric);
//     gmsh::model::setCurrent(metricOptions.modelForMesh);

//     gmsh::model::getEntities(dimTags, 2);
//     if(dimTags.size() > 1) {
//       feWarning("Gmsh model has more than one surface");
//     } else {
//       faceTag = dimTags[0].second;
//     }

//     std::vector<int> viewTags;
//     gmsh::view::getTags(viewTags);
//     feInfo("There are %d views in the gmsh model : ", viewTags.size());
//     for(auto val : viewTags) {
//       feInfo("View %d with index %d", val, gmsh::view::getIndex(val));
//     }

//     metric->setMetricViewTag(viewTags[0]);

//     // computePointsUsingScaledCrossFieldPlanarP2(
//     //   metricOptions.modelForMetric.c_str(), metricOptions.modelForMesh.c_str(),
//     //   metric->getMetricViewTag(), faceTag, pts, errorSquaredCallback, metricOptions.insideCallback,
//     //   nullptr, onlyGenerateVertices, evaluateFieldFromRecoveryCallback, (void *)recovery,
//     // interpolateMetricP1WithDerivativesWrapper, interpolateMetricP1Wrapper,
//     // interpolateMetricAndDerivativeOnP2EdgeWrapper, interpolateMetricP1Wrapper1D,
//     // interpolateMetricAndDerivativeOnP2EdgeWrapper1D,
//     // (void *)metric);

//     gmsh::write(metricOptions.adaptedMeshName);
//   }

// #else
//   printf("In feAdaptMesh : Error - Gmsh is required to generate curved meshes.\n");
// #endif
// }