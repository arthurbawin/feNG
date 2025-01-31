#ifndef _FE_MINIMIZE_INTERPOLATION_ERROR_
#define _FE_MINIMIZE_INTERPOLATION_ERROR_

#include "feMesh.h"
#include "feNorm.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#if defined(GMSH_WITH_CURVED_MESHING)
#include "curvedMesh_structs.h"
#endif
#endif

extern feNewRecovery *activeRecovery;
extern feSpace *activeIntSpace;
extern feNumber *activeNumbering;
extern feSolution *activeSolution;
extern feFunction *activeExactSolution;
extern feVectorFunction *activeExactSolutionGradient; 

extern feMesh *activeMesh;
extern feCncGeo *activeConnectivity;
extern feNorm *activeNorm;

#if defined(HAVE_GMSH)
void computeInterpolationErrorOnEachElement();
double computeInterpolationError(const edgeAndVertexData &data, double *modifiedCoord);
void computeInterpolationErrorGradient(const edgeAndVertexData &data,
                                       const double *modifiedCoord,
                                       double gradient[2]);
double computeInterpolationErrorCallback_EdgesAndVertices(const std::vector<double> &newPositions,
                                                          const curvedAdaptData &adaptData,
                                                          const std::vector<PolyMesh::Vertex*> &verticesToMove,
                                                          const std::vector<PolyMesh::HalfEdge*> &edgesToMove);
void applyCurvatureToFeMesh(const edgeAndVertexData &data, const double *modifiedCoord);
void getMidnodeTags(const SPoint2 edge[2], const double tol, int &elementTag, int &localTag, int &globalTag);
void getPolyMeshVertexTags(const SPoint2 &p, const double tol, 
  std::vector<int> &elementTags, std::vector<int> &localTags, int &globalTag);
#endif

#endif
