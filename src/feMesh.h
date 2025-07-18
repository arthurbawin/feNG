#ifndef _FEMESH_
#define _FEMESH_

#include "feNG.h"
#include "feMessage.h"
#include "feCncGeo.h"
#include "feSpace.h"

#include "feVertex.h"
#include "feEdge.h"
#include "feTriangle.h"
#include "feTetrahedron.h"

#include "rtree.h"

#if defined(HAVE_GMSH) && defined(GMSH_WITH_CURVED_MESHING)
#include "curvedMesh_structs.h"
#endif
//
// Abstract mesh class. Meshes are defined through the derived
// classes according to their highest space dimension:
//  0D mesh: feMesh0DP0
//  1D mesh: feMesh1DP1
//  2D mesh: feMesh2DP1
//
class feMesh
{
protected:
  // Name of the mesh (currently unused)
  std::string _ID;
  // Maximum dimension of the mesh elements
  int _dim;
  // Number of vertices
  int _nVertices;
  // Number of vertices defined on the same geometric entity as another node
  int _nVerticesWithNoPhysical;
  // Number of edges (1D-facets)
  int _nEdges;
  // Number of triangle facets (2D-facets)
  int _nFaces;
  // Total number of elements of all types and all dimensions
  int _nTotalElm;
  // Number of non-boundary elements
  int _nInteriorElm;
  // Number of boundary elements
  int _nBoundaryElm;

  // Highest order of the mesh elements
  int _order;

  // Mesh vertices (all)
  std::vector<Vertex> _vertices;
  // Mesh vertices for each Physical Entity (geometric connectivity)
  // std::map<const feCncGeo *, std::vector<Vertex> > _verticesPerCnc;

  // Number, vector and map of geometric connectivities (domains) in the mesh
  int _nCncGeo;
  std::vector<feCncGeo*> _cncGeo;
  std::map<std::string, int> _cncGeoMap;

  std::vector<const feCncGeo*> _0DConnectivities;
  std::vector<const feCncGeo*> _1DConnectivities;
  std::vector<const feCncGeo*> _2DConnectivities;
  std::vector<const feCncGeo*> _3DConnectivities;

public:
  // Map from non-sequential Gmsh vertices tag to internal sequential tag
  std::map<int, int> _verticesMap;
  // Set of (unique) edges
  std::set<Edge, EdgeLessThan> _edges;
  // Unordered vector of edges
  std::vector<const Edge *> _edgesVec;
  std::map<int, const Edge *> _edgesMap;
  // The P2 edge vertex for each edge, if curved. Empty otherwise.
  std::map<const Edge*, Vertex*> _edge2midnode;
  // The P2 displacement vector xMid - (x0+x1)/2 for each edge
  std::map<const Edge*, double[2]> _edge2alpha;
  // Mesh elements, hardcoded triangles for now
  std::vector<Triangle*> _elements;

///////////////////////////////////////////////
protected:
  // The elements by dimension.
  // The lower dimensional elements are created as boundaries of the higher-dimensional elements.
  std::set<Edge, EdgeLessThan> _allEdges;
  std::set<Triangle, TriangleLessThan> _allTriangles;
  std::set<Tetrahedron, TetrahedronLessThan> _allTetrahedra;

  // // std::vector<Vertex> _vertices_vec;
  // std::vector<Edge> _edges_vec;
  // std::vector<Triangle> _triangles_vec;
  std::vector<Tetrahedron> _tetrahedra;
///////////////////////////////////////////////

public:
  // The constructor of the base class should not be called directly. Instead,
  // call the constructor of the derived class with matching dimension, e.g. mesh2DP1.
  feMesh(int nNod = 0, int dim = 0, int nCncGeo = 0, std::string ID = "")
    : _ID(ID), _dim(dim), _nVertices(nNod), _nEdges(0), _nCncGeo(nCncGeo){};
  virtual ~feMesh() {}

  std::string getID() { return _ID; }
  int getDim() { return _dim; }
  int getNumVertices() { return _nVertices; }
  int getNumVerticesWithoutPhysical() { return _nVerticesWithNoPhysical; }
  int getNumEdges() { return _nEdges; }
  int getNumFaces() { return _nFaces; }
  int getNumElements() { return _nTotalElm; }
  int getNumInteriorElements() { return _nInteriorElm; }
  int getNumBoundaryElements() { return _nBoundaryElm; }
  int getOrder() const { return _order; }
  int getVertexSequentialTagFromGmshTag(int gmshNodeTag) { return _verticesMap[gmshNodeTag]; }
  int getNumCncGeo() { return _nCncGeo; }

  std::vector<feCncGeo*> getCncGeo()
  { 
    std::vector<feCncGeo*> connectivities;
    for(auto const &cnc : _cncGeo)
      connectivities.push_back(cnc);
    return connectivities;
  }
  std::vector<const feCncGeo*> getCncGeo() const
  { 
    std::vector<const feCncGeo*> connectivities;
    for(auto const &cnc : _cncGeo)
      connectivities.push_back(cnc);
    return connectivities;
  }

  std::vector<Vertex> &getVertices() { return _vertices; }
  Vertex *getVertex(int iVertex) { return &_vertices[iVertex]; }
  Vertex *getVertexFromGmshNodeTag(int gmshNodeTag)
  {
    return &_vertices[_verticesMap[gmshNodeTag]];
  }
  const Edge *getEdge(int edgeTag) { return _edgesMap[edgeTag]; }
  const std::map<int, const Edge *> &getEdgesMap() const { return _edgesMap; }
  std::vector<const Edge *> &getEdges() { return _edgesVec; }

  void recomputeDisplacementVectors();
  double getMaxNormEdgeDisplacement();
  double getLpNormEdgeDisplacement(int p);

  // Write in geoCoord the physical coordinates of the nodes on the
  // element with LOCAL tag numElem on the connectivity named
  // cncGeoID or numbered cncGeoTag. The geoCoord vector should have size
  // equal to 3 * numNodesPerElement.
  //
  // The coordinates are structured as follows:
  //   [x1 y1 z1 x2 y2 z2 ... xn yn zn]
  void getCoord(const feCncGeo *cnc, const int numElem, std::vector<double> &geoCoord);
  void getCoord(std::string const &cncGeoID, const int numElm, std::vector<double> &geoCoord);
  void getCoord(const int cncGeoTag, const int numElm, std::vector<double> &geoCoord);

  // Get/set the coordinates of mesh vertex numbered iVertex to coord (size = 3)
  void getVertexCoord(const int iVertex, double coord[3]);
  // Return an error status if the new vertex position causes one or more jacobian
  // determinant to become negative
  feStatus setVertexCoord(const int iVertex, const double coord[3]);

  // Return the tag (unique >=0 integer identifier) of connectivity named cncGeoID
  int getCncGeoTag(std::string const &cncGeoID) const;

  // Return the pointer to the connectivity
  feCncGeo *getCncGeoByName(const std::string &cncGeoID);
  feCncGeo *getCncGeoByTag(const int cncGeoTag);
  const feCncGeo *getCncGeoByName(const std::string &cncGeoID) const;
  const feCncGeo *getCncGeoByTag(const int cncGeoTag) const;

public:
  // Gmsh inspired iterators
  typedef std::vector<const feCncGeo*>::iterator cncIter;
  typedef std::vector<const feCncGeo*>::const_iterator cncConstIter;

public:
  cncIter first0DCnc() { return _0DConnectivities.begin(); }
  cncIter first1DCnc() { return _1DConnectivities.begin(); }
  cncIter first2DCnc() { return _2DConnectivities.begin(); }
  cncIter first3DCnc() { return _3DConnectivities.begin(); }
  cncIter last0DCnc() { return _0DConnectivities.end(); }
  cncIter last1DCnc() { return _1DConnectivities.end(); }
  cncIter last2DCnc() { return _2DConnectivities.end(); }
  cncIter last3DCnc() { return _3DConnectivities.end(); }
  cncConstIter first0DCnc() const { return _0DConnectivities.begin(); }
  cncConstIter first1DCnc() const { return _1DConnectivities.begin(); }
  cncConstIter first2DCnc() const { return _2DConnectivities.begin(); }
  cncConstIter first3DCnc() const { return _3DConnectivities.begin(); }
  cncConstIter last0DCnc() const { return _0DConnectivities.end(); }
  cncConstIter last1DCnc() const { return _1DConnectivities.end(); }
  cncConstIter last2DCnc() const { return _2DConnectivities.end(); }
  cncConstIter last3DCnc() const { return _3DConnectivities.end(); }

  // Return the number of elements on the connectivity 
  int getNumElements(std::string const &cncGeoID);
  int getNumElements(const int cncGeoTag);

  // Return the number of nodes on the reference element
  // associated to the discretization of the connectivity,
  // e.g. 2 for linear edges, 3 for quadratic edges, etc.
  // This applies to the elements used for the geometric
  // discretization, not the interpolation functions (which
  // are in the feSpace).
  int getNumVerticesPerElem(std::string const &cncGeoID);
  int getNumVerticesPerElem(const int cncGeoTag);
  int getNumVerticesPerElemForDimension(const int dim) const;

  // Return the global (unique) tag of element with local
  // tag numElem in matching connectivity.
  int getElement(std::string const &cncGeoID, const int numElem);
  int getElement(const int cncGeoTag, const int numElem);

  // Return the global (unique) tag of the numVertex-th vertex
  // of element with local tag numElem in matching connectivity.
  int getVertex(std::string const &cncGeoID, const int numElem, const int numVertex);
  int getVertex(const int cncGeoTag, const int numElem, const int numVertex);

  // Return the global (unique) tag of the numEdge-th edge
  // of element with local tag numElem in matching connectivity
  int getEdge(std::string const &cncGeoID, const int numElem, const int numEdge);
  int getEdge(const int cncGeoTag, const int numElem, const int numEdge);

  // Return the global (unique) tag of the numFace-th facet (triangle)
  // of element with local tag numElem in matching connectivity
  int getFace(std::string const &cncGeoID, const int numElem, const int numFace);
  int getFace(const int cncGeoTag, const int numElem, const int numFace);

  // Return a pointer to the FE space of matching connectivity
  feSpace *getGeometricSpace(std::string const &cncGeoID);
  feSpace *getGeometricSpace(int cncGeoTag);

  const std::vector<Triangle*> &getTriangles() { return _elements; };

  // Locate the physical node x in the mesh. Return true if the
  // node was found, false otherwise. If the node was found,
  // the parameters are modified as follow:
  //
  //    iElm: the global tag of the element in which the node was found
  //    u: the rst (or xi-eta-zeta) coordinates of x on the reference element
  //      with u[0] = r, u[1] = s, u[2] = t
  //    tol: tolerance (applies only to 2D+ search in tree)
  //
  // If returnLocalElmTag = true, iElm is the local element number in the
  // connectivity on which the vertex was found.
  //
  // If targetConnectivity is the name of an existing 2D connectivity in
  // the mesh, we look for the vertex on this connectivity only, and 
  // return false if it was not found on this one (even if it can be found on
  // another connectivity).
  virtual bool locateVertex(const double *x, int &iElm, double *u, double tol = 1e-12,
                            bool returnLocalElmTag = false, std::string targetConnectivity = "") const = 0;
  virtual bool locateVertexInElements(const feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol = 1e-12) = 0;
};

//
// Trivial 0-dimensional discretization of an point.
// Used to solve ODEs.
//
class feMesh0DP0 : public feMesh
{
protected:
  int _nElm;
  double _xA;
  std::string _domID;
  int _nElmDomain;
  int _nElmBoundary;
  int _nNodDomain;
  int _nNodBoundary;

public:
  feMesh0DP0(double xA, int nElm, std::string domID);
  virtual ~feMesh0DP0();

  bool locateVertex(const double *x, int &iElm, double *u, double tol = 1e-12,
                            bool returnLocalElmTag = false, std::string targetConnectivity = "") const
  {
    UNUSED(x,tol,returnLocalElmTag,targetConnectivity);
    iElm = 0;
    u[0] = 1.;
    return true;
  };
  bool locateVertexInElements(const feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol = 1e-12)
  {
    UNUSED(cnc,x,elementsToSearch,tol);
    iElm = 0;
    u[0] = 1.;
    return true;
  }
};

//
// 1-dimensional uniform discretization of an interval [xBegin, xEnd]
//
class feMesh1DP1 : public feMesh
{
protected:
  // Beginning and end of the interval
  double _xA, _xB;
  // Name of the boundaries and of the interior domain
  std::string _bndA_ID, _bndB_ID, _domID;
  //
  int _nElmPerBoundary;
  int _nNodDomain;
  int _nNodBoundary;

public:
  // Create a 1D mesh of the interval [xBegin,xEnd], with numElements
  // uniform line elements and two boundary point elements.
  // The boundaries are named boundaryBeginName and boundaryEndName,
  // and the interior domain is named domainName.
  feMesh1DP1(const double xBegin, const double xEnd, const int numElements,
             const std::string &boundaryBeginName, const std::string &boundaryEndName,
             const std::string &domainName);
  virtual ~feMesh1DP1();

  bool locateVertex(const double *x, int &iElm, double *u, double tol = 1e-12,
                    bool returnLocalElmTag = false, std::string targetConnectivity = "") const;
  bool locateVertexInElements(const feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol = 1e-12);
};

class feMetaNumber;
class feSolutionContainer;
class feMetric;
class feMetricOptions;
class feRecovery;
class feNewRecovery;

typedef struct rtreeSearchCtxStruct {
  const std::vector<Triangle*> *elements;
  double uvw[3];
  double r[3];
  double x[3];
  double min[2];
  double max[2];
  double tolerance;
  int iElm;
  int iElmLocal;
  int physicalTag;
  bool isFound;
  bool enforceConnectivity = false;
  int targetPhysicalTag;
  int numFound;
} rtreeSearchCtx;

#if !defined(GMSH_WITH_CURVED_MESHING)
enum class curveToMinimize
{ 
  LENGTH, INTERPOLATION_ERROR, DISTORTION
};
#endif

//
// 2-dimensional mesh
//
class feMesh2DP1 : public feMesh
{
private:
  // A Geometric Entity.
  // These are Gmsh's elementary geometric blocks, such as Points,
  // Lines, Surfaces, etc. that are grouped together to form
  // Physical Entities (named domains of subdomains).
  // Each Geometric Entity in the mesh must be part of one (and
  // only one) Physical Enttity (subdomains cannot overlap).
  typedef struct entityStruct {
    int dim;
    int tag;
    int tagGmsh;
    int numPhysicalTags;
    std::vector<int> physicalTags;
    int gmshType; // Element type in Gmsh
    geometricInterpolant interp = geometricInterpolant::NONE;
    int nElm;
    int nNodePerElem = 0;
    int nEdgePerElem;
    int nTriPerElem;
    int nTetPerElem;
    // The global tag of the lower-dimensional boundary elements
    // for the elements on this entity
    std::vector<int> connecNodes;
    std::vector<int> connecEdges;
    std::vector<int> connecTri;
    // The global tag of the elements on this entity,
    // regardless of their dimension.
    std::vector<int> connecElem;
    std::vector<int> isBoundedBy;
    std::vector<int> isBoundaryOf;

    std::vector<int> connecEdgesDebug;
    std::vector<int> connecTrianglesDebug;
  } entity;

  // FIXME: Should probably simply be a feCncGeo, which
  // contains the same (and more) information on a physical entity.
  //
  // A Physical Entity.
  // A group of Geometric Entities of same dimension, with
  // a name. Physical Entities cannot overlap.
  // Equations and boundary conditions are defined on Physical
  // Entities.
  typedef struct physicalEntityStruct {
    std::string name;
    int dim;
    int tag;
    std::vector<int> listEntities; // Tag des entities
    geometryType geometry;
    geometricInterpolant interp = geometricInterpolant::NONE;
    feSpace *geoSpace;
    int nElm;
    int nNodePerElem;
    int nEdgePerElem;
    int nTriPerElem;
    int nTetPerElem;
    std::vector<int> connecNodes;
    std::vector<int> connecEdges;
    std::vector<int> connecTri;
    std::vector<int> connecElem;
    std::vector<int> isBoundedBy;
    std::vector<int> isBoundaryOf;
  } physicalEntity;

  // RTrees and context to quickly search for nodes inside the mesh
  mutable RTree<int, double, 3> _rtree;
  mutable RTree<int, double, 2> _rtree2d;

public:
  // FIXME : Choose more explicit name, for example "physicalEntitiesDescription"
  typedef std::map<std::pair<int, int>, std::string> mapType;

protected:
  // Geometric Entities sorted by pairs <dimension, tag>
  std::map<std::pair<int, int>, entity> _entities;

  // Number of Physical Entities in the mesh
  int _nPhysicalEntities;

  // Physical Entities sorted by pairs <dimension, tag>
  std::map<std::pair<int, int>, physicalEntity> _physicalEntities;

  // A somewhat redundant description of the Physical Entities:
  // <dimension, tag> : name of the Physical Entity
  // Used during mesh adaptation, where the Physical Entities
  // might be lost after remeshing.
  mapType _physicalEntitiesDescription;

  std::vector<int> _sequentialNodeToGmshNode;

  // The number of Gmsh's Points, Curves, Surfaces and Volumes
  // Currently unused outside of feMeshRead.
  int _nPoints;
  int _nCurves;
  int _nSurfaces;
  int _nVolumes;

  // The total number of element (redundant member)
  int _nElm;

  // Informations from the MSH mesh file
  bool _isBinary;
  double _gmshVersion;

public:
  // Empty constructor: the mesh must be created afterwards with feCheck(readGmsh(...))
  feMesh2DP1(){};

  // Create a 2D mesh from GMSH file "meshName".
  // Only reads .msh files with version 4+.
  //
  //  meshName: the relative or absolute path to the mesh file,
  //            e.g. "square.msh" or "../data/square.msh".
  //  curved: true if using a high-order mesh (currently supporting
  //          only 1st and 2nd order geometries), false otherwise.
  //  reversed: true if the mesh should be reversed (vertices are
  //            numbered clockwise instead of counterclockwise).
  //  physicalEntitiesDescription: a map giving the name of the
  //                               Physical Entities for each pair
  //                               <dimension, tag>.
  feMesh2DP1(const std::string &meshName, const bool curved = false, const bool reversed = false,
             const mapType &physicalEntitiesDescription = mapType());
  ~feMesh2DP1();

private:
  // Read version 2.2 or 4+ Gmsh mesh file. Only called in constructor.
  feStatus readMsh2(std::istream &input, const bool curved, const bool reversed,
                    mapType physicalEntitiesDescription);
  feStatus readMsh4(std::istream &input, const bool curved, const bool reversed,
                    mapType physicalEntitiesDescription);
  feStatus readGmsh(const std::string meshName, const bool curved = false,
                    const bool reversed = false,
                    const mapType physicalEntitiesDescription = mapType());

public:
  bool isMeshFileBinary() { return _isBinary; }
  double getGmshVersion() { return _gmshVersion; }

  mapType getPhysicalEntitiesDescription() { return _physicalEntitiesDescription; }

  int getGmshNodeTag(int iVertex) { return _sequentialNodeToGmshNode[iVertex]; }

  double computeAverageSliverness(const double LpNorm);

  bool locateVertex(const double *x, int &iElm, double *u, double tol = 1e-12,
                    bool returnLocalElmTag = false, std::string targetConnectivity = "") const;
  bool locateVertexInElements(const feCncGeo *cnc, const double *x, const std::vector<int> &elementsToSearch, int &iElm,
                                      double *u, double tol = 1e-12);

  // Transfer whole FE solution(s) associated to the current mesh to another mesh.
  // The solution(s) are stored in solutionContainer and are transferred in place.
  feStatus transfer(feMesh2DP1 *otherMesh, feMetaNumber *myMN, feMetaNumber *otherMN,
                    feSolutionContainer *solutionContainer, const std::vector<feSpace *> &mySpaces,
                    const std::vector<feSpace *> &mySpacesEssBC,
                    const std::vector<feSpace *> &otherSpaces);

  feStatus adapt(feMetric *metricField);
  // Default argument is a feMetricOptions initialized by default, see feAdaptMesh.cpp
  feStatus adapt(std::vector<feNewRecovery*> recoveredFields,
                 feMetricOptions &options, 
                 const std::vector<feSpace*> &spaces,
                 const std::vector<feSpace*> &essentialSpaces,
                 feSpace *spaceForAdaptation,
                 feSolution *discreteSolution,
                 feFunction *exactSolution,
                 feVectorFunction *exactGradient,
                 bool curve,
                 bool isBackmeshP2,
                 bool curveMMGmesh,
                 curveToMinimize target,
                 bool generateAnisoMeshBeforeCurving
#if defined(HAVE_GMSH) && defined(GMSH_WITH_CURVED_MESHING)
                 ,
                 directionField targetDirectionField,
                 vertexSpawning targetVertexSpawning
#endif
                 );

  void drawConnectivityToPOSfile(const std::string &cncName, const std::string &fileName);
};

#endif
