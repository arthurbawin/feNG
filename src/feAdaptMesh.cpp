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

  feInfo("Running MMG with: %s", command.data());

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
                           bool curve, bool isBackmeshP2, bool curveMMGmesh,
                           curveToMinimize target,
                           bool generateAnisoMeshBeforeCurving
#if defined(HAVE_GMSH)
                           ,
                           directionField targetDirectionField,
                           vertexSpawning targetVertexSpawning
#endif
                 )
{
  if(curve && !isBackmeshP2) {
    return feErrorMsg(FE_STATUS_ERROR, "Backmesh must be P2 when curving the mesh.");
  }

#if defined(HAVE_GMSH)

  // Step 1: Compute the metric tensor field and write metrics to options.mmgInputMeshfile
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
  gmsh::model::getCurrent(options.modelForMetric);

  feMetric metricField(recoveredFields, options);
  metricField._backmeshOrder = getGeometricInterpolantDegree(recoveredFields[0]->_cnc->getInterpolant());
  metricField._nVerticesPerElmOnBackmesh = recoveredFields[0]->_cnc->getFeSpace()->getNumFunctions();

  tic();
  feStatus s = metricField.computeMetrics();
  options.userValue = metricField._options.userValue;
  feInfoCond(FE_VERBOSE > 0, "\t\tComputed metric tensors in %f s", toc());
  if(s != FE_STATUS_OK) {
    gmsh::finalize();
    gmshWasInitialized = false;
    return s;
  }

  // To test interpolation error only
  // return FE_STATUS_OK;

  // Step 2: Create aniso mesh

  // If we curve to minimize interpolation error, we have to use the feMesh
  // built on the current mesh. So we don't generate a new anisotropic mesh
  // before curving.
  if(generateAnisoMeshBeforeCurving && !(curve && target == curveToMinimize::INTERPOLATION_ERROR))
  {
    // Save directly to .msh to preserve Physical Entities
    std::string cmd1 = "mmg2d_O3 " + options.mmgInputMeshfile + " -v 10 -hgrad -1 -o " + options.mmgOutputMeshfile;
    if(FE_VERBOSE == VERBOSE_NONE) {
      // Write MMG console outputs to logMMG.txt
      cmd1 += " > logMMG.txt";
    }

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

    // if(metricField._backmeshOrder > 1) {
    //   gmsh::model::mesh::setOrder(metricField._backmeshOrder);
    // }

    // Step 3: Check physical entities
    gmsh::vectorpair physicalGroups;
    gmsh::vectorpair allGeometricEntities;
    gmsh::model::getPhysicalGroups(physicalGroups);
    gmsh::model::getEntities(allGeometricEntities);

    std::map<std::pair<int, int>, std::vector<int>> entitiesForPhysical;

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

    // Petite astuce: écrire au format 2.2 pour éliminer des
    // entités géométriques douteuses 0D de MMG...
    gmsh::write(options.adaptedMeshName);
    gmsh::clear();
    gmsh::open(options.adaptedMeshName);

    // gmsh::model::mesh::reverse();

    // Write P2 aniso mesh
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::model::mesh::setOrder(metricField._backmeshOrder);
    gmsh::write("anisoadapted.msh");
    gmsh::write(options.adaptedMeshName);
    gmsh::model::mesh::setOrder(1);
  }

  if(curve) {

    // Must merge the local Gmsh files first
    // Quick fix in the meantime
    #if defined(GMSH_WITH_CURVED_MESHING)
    
    feNorm *norm;
    if(target == curveToMinimize::INTERPOLATION_ERROR) {
      // Tools to compute interpolation error while curving
      createNorm(norm, L2_ERROR, {spaceForAdaptation}, discreteSolution, exactSolution, exactGradient);
      norm->setRecovery(recoveredFields[0]); 
      activeMesh = this;
      activeRecovery = recoveredFields[0];
      activeNorm = norm;
      activeConnectivity = recoveredFields[0]->_cnc;
      activeIntSpace = spaceForAdaptation;
      activeSolution = discreteSolution;
      activeExactSolution = exactSolution;
      // activeExactSolution = nullptr;
      activeExactSolutionGradient = exactGradient;
      // activeExactSolutionGradient = nullptr;
      computeInterpolationErrorOnEachElement();
    }

    if(metricField._backmeshOrder > 1) {
      // GFace2PolyMesh only takes P1 meshes
      gmsh::model::mesh::setOrder(1);
    }

    // For now only curve a single 2D surface
    int faceTag;
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags, 2);
    if(dimTags.size() > 1) {
      for(auto pair : dimTags) {
        feInfo("Gmsh 2D entity: dim = %d - tag = %d", pair.first, pair.second);
      }
      return feErrorMsg(FE_STATUS_ERROR, "Gmsh model has more than one surface");
    } else {
      faceTag = dimTags[0].second;
    }

    // std::vector<int> viewTags;
    // gmsh::view::getTags(viewTags);
    // metricField.setMetricViewTag(viewTags[0]);

    // Check there is an "inside" callback
    if(options.insideCallback == nullptr){
      return feErrorMsg(FE_STATUS_ERROR,
        "Provide an \"inside\" callback to peform curved mesh adaptation");
    }

    curvedMeshOptions meshOptions;

    std::string modelName;
    gmsh::model::getCurrent(modelName);

    // modelForMetric, modelForMesh and VIEW_TAG are deprecated
    // meshOptions.modelForMetric = options.modelForMetric.data();
    meshOptions.modelForMetric = modelName.data();
    // meshOptions.modelForMesh = options.modelForMetric.data();
    meshOptions.modelForMesh = modelName.data();
    meshOptions.VIEW_TAG = metricField.getMetricViewTag();

    meshOptions.targetDirectionField = targetDirectionField;
    meshOptions.targetVertexSpawning = targetVertexSpawning;

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
    meshOptions.interpolateMetricP2CallbackWithoutDerivatives = interpolateMetricP2CallbackWithoutDerivatives;
    meshOptions.interpolateMetricP2CallbackWithoutDerivativesExplicit = interpolateMetricP2CallbackWithoutDerivativesExplicit;
    meshOptions.interpolateMetricP2CallbackWithDerivatives = interpolateMetricP2CallbackWithDerivatives;
    meshOptions.interpolateMetricP2CallbackLog = interpolateMetricP2CallbackLog;

    tic();
    createCurvedMesh(meshOptions);
    feInfoCond(FE_VERBOSE > 0, "\t\tGenerated curved mesh in %f s", toc());

    // Metric interpolation error
    // Or
    // Final error after curving when minimizing interpolation error
    options.userValue = meshOptions.userValue;

    // gmsh::model::setCurrent(options.modelForMetric.data());
    gmsh::option::setNumber("Mesh.MshFileVersion", 4.1);
    gmsh::write(options.adaptedMeshName);
    // system("gmsh adapted.msh");
    #else
      return feErrorMsg(FE_STATUS_ERROR,
                          "Cannot generate curved mesh with this version of Gmsh."
                          "Merge branch first.");
    #endif
  }

  gmsh::clear();
  gmsh::finalize();

  gmshWasInitialized = false;
#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to compute metric tensors and adapt the mesh!");
#endif
  return FE_STATUS_OK;
}

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