#include "feNG.h"
#include "feMetric.h"
#include "feMetricTools.h"
#include "feNewRecovery.h"
#include "feMatrixInterface.h"

extern int FE_VERBOSE;

#define EXPORT_FIELDS_TO_GMSH

#define XLOC_TOL 1e-16
#define ORIGINAL_ZHANG_NAGA
#define APPLY_SCALING
// #define TREAT_BOUNDARIES
// #define AVERAGE_EVALUATIONS_FROM_DOFS

// The names of the reconstructed fields
static std::map<std::pair<int, int>, std::string> suffix = {
  {{0, 0}, ""},

  {{1, 0}, "dx"},     {{1, 1}, "dy"},

  {{2, 0}, "dxx"},    {{2, 1}, "dxy"},
  {{2, 2}, "dyx"},    {{2, 3}, "dyy"},

  {{3, 0}, "dxxx"},   {{3, 1}, "dxxy"},
  {{3, 2}, "dxyx"},   {{3, 3}, "dxyy"},
  {{3, 4}, "dyxx"},   {{3, 5}, "dyxy"},
  {{3, 6}, "dyyx"},   {{3, 7}, "dyyy"},

  {{4, 0}, "dxxxx"},  {{4, 1}, "dxxxy"},
  {{4, 2}, "dxxyx"},  {{4, 3}, "dxxyy"},
  {{4, 4}, "dxyxx"},  {{4, 5}, "dxyxy"},
  {{4, 6}, "dxyyx"},  {{4, 7}, "dxyyy"},
  {{4, 8}, "dyxxx"},  {{4, 9}, "dyxxy"},
  {{4, 10}, "dyxyx"}, {{4, 11}, "dyxyy"},
  {{4, 12}, "dyyxx"}, {{4, 13}, "dyyxy"},
  {{4, 14}, "dyyyx"}, {{4, 15}, "dyyyy"},

  {{5, 0}, "dxxxxx"},  {{5, 1}, "dxxxxy"},
  {{5, 2}, "dxxxyx"},  {{5, 3}, "dxxxyy"},
  {{5, 4}, "dxxyxx"},  {{5, 5}, "dxxyxy"},
  {{5, 6}, "dxxyyx"},  {{5, 7}, "dxxyyy"},
  {{5, 8}, "dxyxxx"},  {{5, 9}, "dxyxxy"},
  {{5, 10}, "dxyxyx"}, {{5, 11}, "dxyxyy"},
  {{5, 12}, "dxyyxx"}, {{5, 13}, "dxyyxy"},
  {{5, 14}, "dxyyyx"}, {{5, 15}, "dxyyyy"},
  {{5, 16}, "dyxxxx"}, {{5, 17}, "dyxxxy"},
  {{5, 18}, "dyxxyx"}, {{5, 19}, "dyxxyy"},
  {{5, 20}, "dyxyxx"}, {{5, 21}, "dyxyxy"},
  {{5, 22}, "dyxyyx"}, {{5, 23}, "dyxyyy"},
  {{5, 24}, "dyyxxx"}, {{5, 25}, "dyyxxy"},
  {{5, 26}, "dyyxyx"}, {{5, 27}, "dyyxyy"},
  {{5, 28}, "dyyyxx"}, {{5, 29}, "dyyyxy"},
  {{5, 30}, "dyyyyx"}, {{5, 31}, "dyyyyy"},
};

static inline double myPow(double base, int exp)
{
  // return pow(base, exp);
  double res = 1.;
  for(int i = 0; i < exp; ++i) res *= base;
  return res;
}

static void drawCircle(FILE *f, double *x, double radius, int nPoints)
{
  int nt = nPoints;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  double dtheta = (2. * M_PI) / nPoints;
  for(int i = 0; i < nt; ++i) {
    double theta = i * dtheta;
    xP[i] = x[0] + radius * cos(theta);
    yP[i] = x[1] + radius * sin(theta);
  }

  for(int i = 0; i < nt; ++i) {
    if(i != nt - 1) {
      fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[i + 1],
              yP[i + 1], 0., 1, 1);
    } else {
      fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[0],
              yP[0], 0., 1, 1);
    }
  }
}

feNewPatch::feNewPatch(const feCncGeo *cnc, feMesh *mesh, bool reconstructAtHighOrderNodes,
                       bool useOriginalZhangNagaPatchDefinition, int degreeSolution,
                       feMetric *metricField)
{
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getVerticesConnectivity();

  _degSol = degreeSolution;

  int nElm = cnc->getNumElements();
  _nNodePerElm = cnc->getNumVerticesPerElem();
  _nVertices = _vertices.size() / _nNodePerElm;
  _nEdgePerElm = cnc->getNumEdgesPerElem();

  // If mesh if high-order, keep only the P1 vertices.
  // This is because the Zhang-Naga method reconstructs fields at vertices,
  // then extends to edges and inside by evaluating/averaging.
  _elementOrder = getGeometricInterpolantDegree(cnc->getInterpolant());

  _reconstructAtHighOrderNodes = reconstructAtHighOrderNodes;

  // Do not fit polynomials at high order nodes if interpolant is linear
  // because there is no information at edge mesh vertices
  if(_degSol == 1)
    _reconstructAtHighOrderNodes = false;

  // We have to reconstruct at P2 nodes to compute metrics on
  // the P2 backmesh
  std::map<int, std::pair<int,int>> edgeVertexToAdjacentVertices;
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < 3; ++j) {
      isOnEdge.insert({_vertices[_nNodePerElm * i + j], false});
    }
    for(int j = 3; j < _nNodePerElm; ++j) {
      isOnEdge.insert({_vertices[_nNodePerElm * i + j], true});
      int v1 = _vertices[_nNodePerElm * i + j - 3];
      int v2 = _vertices[_nNodePerElm * i + (j - 2) % 3];
      edgeVertexToAdjacentVertices[_vertices[_nNodePerElm * i + j]] = {v1,v2};
    }
  }

  if(!reconstructAtHighOrderNodes) {
    // Vertices to remove are replaced by the first vertex tag,
    // since the vector will be sorted and uniqued.
    // Only keep the first 3 vertices.
    int firstVertex = _vertices[0];
    for(int i = 0; i < nElm; ++i) {
      for(int j = 3; j < _nNodePerElm; ++j) {
        _vertices[_nNodePerElm * i + j] = firstVertex;
      }
    }
  }
 
  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  const std::vector<int> &connecNodes = cnc->getVerticesConnectivity();

  // Patches of adjacent elements
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < _nNodePerElm; ++j) {
      vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i);
      vertToElems_singleLayer[connecNodes[_nNodePerElm * i + j]].insert(i);
      _isVertexBoundary[connecNodes[_nNodePerElm * i + j]] = false;
    }
  }

  for(auto p : vertToElems_singleLayer) {
    size_t v = p.first;
    for(auto elm : p.second) {
      for(int j = 0; j < _nNodePerElm; ++j) {
        vertToVerts_singleLayer[v].insert(connecNodes[_nNodePerElm * elm + j]);
      }
    }
  }

  // Mark boundary vertices
  // Vertices on cnc with dim = dimMesh - 1 are considered boundary
  for(auto cnc : mesh->getCncGeo()) {
    if(cnc->getDim() == mesh->getDim() - 1) {
      for(auto bnd_vertex : cnc->getVerticesConnectivity()) {
        _isVertexBoundary.at(bnd_vertex) = true;
      }
    }
  }

  std::vector<double> geoCoord(3 * _nNodePerElm, 0.);

  // Dimension of the polynomial basis for the recovery
  int basisDimension = (_degSol + 2) * (_degSol + 3) / 2;

  // int minVertices = basisDimension;
  // int minVertices = round(1.5 * basisDimension);
  int minVertices = 2 * basisDimension;

  if(_nVertices < minVertices) {
    feErrorMsg(FE_STATUS_ERROR,
               "Mesh has too few vertices (%d) to compute solution recovery. Need at least %d "
               "vertices to compute recovery for interior nodes.",
               _nVertices, minVertices);
    exit(-1);
  }

  _useOriginalZhangNagaPatchDefinition = useOriginalZhangNagaPatchDefinition;

  if(_useOriginalZhangNagaPatchDefinition) {

    int minVerticesOnBoundary = minVertices;

    // Compute length of incident edges for original patch
    for(int i = 0; i < nElm; ++i) {
      for(int j = 0; j < 3; ++j) {
        _hMax[connecNodes[_nNodePerElm * i + j]] = 0.;
      }
    }

    for(int i = 0; i < nElm; ++i) {
      mesh->getCoord(cnc->getTag(), i, geoCoord);
      double x0 = geoCoord[0];
      double y0 = geoCoord[1];
      double x1 = geoCoord[3];
      double y1 = geoCoord[4];
      double x2 = geoCoord[6];
      double y2 = geoCoord[7];
      // Use Euclidian metric (bad for aniso meshes?)
      double l0 = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
      double l1 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
      double l2 = sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0));
      double L[3] = {l0, l1, l2};

      // For P1 vertices take the longest edge
      for(int j = 0; j < 3; ++j) {
        int vertex = connecNodes[_nNodePerElm * i + j];
        double longestIncident = fmax(L[j], L[(j + 2) % 3]);
        _hMax.at(vertex) = fmax(_hMax.at(vertex), longestIncident);
      }

      // For P2 midnodes, take the longest (straight) vertex-to-vertex distance in the triangle
      // for(int j = 3; j < _nNodePerElm; ++j) {
      //   int vertex = connecNodes[_nNodePerElm * i + j];
      //   double xj = geoCoord[3 * j + 0];
      //   double yj = geoCoord[3 * j + 1];
      //   l0 = sqrt((x0 - xj) * (x0 - xj) + (y0 - yj) * (y0 - yj));
      //   l1 = sqrt((x1 - xj) * (x1 - xj) + (y1 - yj) * (y1 - yj));
      //   l2 = sqrt((x2 - xj) * (x2 - xj) + (y2 - yj) * (y2 - yj));
      //   double longest = fmax(l0, fmax(l1, l2));
      //   _hMax.at(vertex) = fmax(_hMax.at(vertex), longest);
      // }
    }

    // Or compute longest connected edge at midnodes as the max of the hMax of extremities?
    // for(int i = 0; i < nElm; ++i) {
    //   for(int j = 3; j < 6; ++j) {
    //     int e1 = j - 3;
    //     int e2 = (j - 2) % 3;
    //     _hMax.at(connecNodes[_nNodePerElm * i + j]) = fmax(
    //       _hMax.at(connecNodes[_nNodePerElm * i + e1]), _hMax.at(connecNodes[_nNodePerElm * i +
    //       e2]));
    //   }
    // }

    // Debug check: plot circle of radius hMax around each vertex
    if(false) {
      for(auto &p : _hMax) {
        std::string isIso = metricField ? "aniso" : "iso";
        std::string posfile = "cercle" + std::to_string(p.first) + isIso + ".pos";
        FILE *file = fopen(posfile.data(), "w");
        fprintf(file, "View\"cercle%s%d\"{\n", isIso.data(), p.first);
        double x[2] = {0., 0.};
        x[0] = mesh->getVertex(p.first)->x();
        x[1] = mesh->getVertex(p.first)->y();
        if(metricField) {
          drawSingleEllipse(file, x, metricField->getMetricAtSequentialTag(p.first), 1, 30);
        } else {
          double r = p.second;
          drawCircle(file, x, r, 30);
        }
        fprintf(file, "};\n");
        fclose(file);
      }
    }

    // Dumb brute-force: check all vertices and add to patch if closer than _hMax
    for(auto &pair : _hMax) {
      int vertex = pair.first;
      Vertex *vi = mesh->getVertex(vertex);

      bool enoughAdjacentVertices = false;

      do {
        for(auto vertex_j : _vertices) {
          Vertex *vj = mesh->getVertex(vertex_j);
          double dist = sqrt((vi->x() - vj->x()) * (vi->x() - vj->x()) +
                             (vi->y() - vj->y()) * (vi->y() - vj->y()));
          if(dist <= 1.01 * pair.second) {
            vertToVerts[vertex].insert(vertex_j);
          }
        }

        pair.second *= 2.;

        if(_isVertexBoundary.at(vertex) && vertToVerts.at(vertex).size() >= minVerticesOnBoundary)
          enoughAdjacentVertices = true;
        if(!_isVertexBoundary.at(vertex) && vertToVerts.at(vertex).size() >= minVertices)
          enoughAdjacentVertices = true;

      } while(!enoughAdjacentVertices);
    }

    if(_reconstructAtHighOrderNodes) {
      // Another strategy for P2 patches: just merge the patches of extremities
      for(int i = 0; i < nElm; ++i) {
        for(int j = 3; j < 6; ++j) {
          int e1 = j - 3;
          int e2 = (j - 2) % 3;
          int vertex = connecNodes[_nNodePerElm * i + j];
          int v1 = connecNodes[_nNodePerElm * i + e1];
          int v2 = connecNodes[_nNodePerElm * i + e2];
          for(auto connectedVertex : vertToVerts.at(v1))
            vertToVerts[vertex].insert(connectedVertex);
          for(auto connectedVertex : vertToVerts.at(v2))
            vertToVerts[vertex].insert(connectedVertex);
        }
      }
    }

  } else {
    // Do not use patch definition from Zhang & Naga,
    // instead use only layers of elements and associated vertices.
    // Maybe more suited for anisotropic meshes.

    std::set<int> toAdd;

    // Essai en séparant intérieur et frontiere
    // D'abord déterminer les patches interieures
    // Elles doivent contenir au moins "minVertices" sommets
    // On ajoute juste des couches d'éléments pour y arriver

    // Also, check that there is at least one inner P1 vertex in the mesh
    bool OK = false;
    for(auto &pair : vertToElems) {
      int vertex = pair.first;
      if(!_isVertexBoundary.at(vertex) && !isOnEdge.at(vertex)) {

        OK = true;
        bool enoughAdjacentVertices = false;

        do {
          // Insert vertices of adjacent elements to set of connected vertices
          for(auto elem : pair.second) {
            int nVerticesPerElem = _reconstructAtHighOrderNodes ? _nNodePerElm : 3;
            for(int j = 0; j < nVerticesPerElem; ++j) {
              vertToVerts[vertex].insert(connecNodes[_nNodePerElm * elem + j]);
              if(_isVertexBoundary[vertex] &&
                 !_isVertexBoundary[connecNodes[_nNodePerElm * elem + j]])
                vertToInteriorVerts[vertex].insert(connecNodes[_nNodePerElm * elem + j]);
            }
          }

          // Check that there are enough adjacent vertices
          if(!_isVertexBoundary[vertex] && vertToVerts[vertex].size() >= minVertices)
            enoughAdjacentVertices = true;

          if(!enoughAdjacentVertices) {
            toAdd.clear();

            // Add one layer of elements
            // These are the elements connected to the current patch of elements
            // Take only the first layer, as the patch of elements may have grown already
            for(auto currentElem : pair.second) {
              for(int j = 0; j < 3; ++j) {
                for(auto elemToAdd :
                    vertToElems_singleLayer[connecNodes[_nNodePerElm * currentElem + j]]) {
                  toAdd.insert(elemToAdd);
                }
              }
            }

            for(auto elemToAdd : toAdd) {
              pair.second.insert(elemToAdd);
            }
          }
        } while(!enoughAdjacentVertices);
      }
    }

    if(!OK) {
      // All P1 vertices are on the boundary. Need to look into it, if necessary.
      feErrorMsg(FE_STATUS_ERROR,
               "Mesh has no inner P1 vertex. All P1 vertices are on the boundary."
               " Cannot construct patches for derivatives recovery for now.");
      exit(-1);
    }

    // Boundary patches
    int nVerticesPerElem = _reconstructAtHighOrderNodes ? _nNodePerElm : 3;
    for(auto v : _vertices) {
      if(_isVertexBoundary.at(v) && !isOnEdge.at(v)) {

        // Start with the first layer
        boundaryVertToElems[v] = vertToElems_singleLayer.at(v);

        // Check if there is an interior vertex connected
        bool OK = false;
        for(auto connectedElem : boundaryVertToElems.at(v)) {
          for(int i = 0; i < 3; ++i) {
            int connectedVertex = connecNodes[_nNodePerElm * connectedElem + i];
            if(!_isVertexBoundary.at(connectedVertex)) {
              OK = true;
              break;
            }
          }
        }

        while(!OK) {
          // Add a layer of elements
          toAdd.clear();
          for(auto connectedElem : boundaryVertToElems.at(v)) {
            for(int i = 0; i < 3; ++i) {
              int connectedVertex = connecNodes[_nNodePerElm * connectedElem + i];
              for(auto newElem : vertToElems_singleLayer.at(connectedVertex)) {
                toAdd.insert(newElem);
              }
            }
          }
          for(auto elemToAdd : toAdd) {
            boundaryVertToElems.at(v).insert(elemToAdd);
          }

          // Check if there is an interior vertex connected
          for(auto connectedElem : boundaryVertToElems.at(v)) {
            for(int i = 0; i < 3; ++i) {
              int connectedVertex = connecNodes[_nNodePerElm * connectedElem + i];
              if(!_isVertexBoundary.at(connectedVertex)) {
                OK = true;
                break;
              }
            }
          }
        }

        // Add the element patch of each interior vertex in the current vertex patch of the boundary vertex
        toAdd.clear();
        for(auto elem : boundaryVertToElems.at(v)) {
          for(int j = 0; j < nVerticesPerElem; ++j) {
            int connectedVertex = connecNodes[_nNodePerElm * elem + j];
            if(!_isVertexBoundary.at(connectedVertex)) {
              for(auto newElem : vertToElems.at(connectedVertex))
                toAdd.insert(newElem);
            }
          }
        }

        // Add the first layer of the boundary vertices in the current patch
        for(auto elem : boundaryVertToElems.at(v)) {
          for(int j = 0; j < nVerticesPerElem; ++j) {
            int connectedVertex = connecNodes[_nNodePerElm * elem + j];
            if(_isVertexBoundary.at(connectedVertex)) {
              for(auto newElem : vertToElems_singleLayer.at(connectedVertex))
                toAdd.insert(newElem);
            }
          }
        }
        for(auto elemToAdd : toAdd) {
          boundaryVertToElems.at(v).insert(elemToAdd);
        }

        // Create the vertToVert map at boundary vertices
        for(auto elem : boundaryVertToElems.at(v)) {
          for(int j = 0; j < nVerticesPerElem; ++j) {
            vertToVerts[v].insert(connecNodes[_nNodePerElm * elem + j]);
          }
        }

        vertToElems.at(v) = boundaryVertToElems.at(v);
      }
    }

    // Finally create edge patches as the union of both extremities
    for(auto v : _vertices) {
      if(isOnEdge.at(v)) {

        // Add the vertex patch of each interior vertex in the current vertex patch of the boundary vertex
        toAdd.clear();
        auto pair = edgeVertexToAdjacentVertices.at(v);
        for(auto newVertex : vertToVerts.at(pair.first))
          toAdd.insert(newVertex);
        for(auto newVertex : vertToVerts.at(pair.second))
          toAdd.insert(newVertex);
        for(auto vertexToAdd : toAdd) {
          vertToVerts[v].insert(vertexToAdd);
        }

        vertToElems.at(v).clear();
        for(auto elem : vertToElems.at(pair.first))
          vertToElems.at(v).insert(elem);
        for(auto elem : vertToElems.at(pair.second))
          vertToElems.at(v).insert(elem);

      }
    }

  }

  // Compute the scaling
  for(auto v : _vertices) {
    Vertex *vi = mesh->getVertex(v);
    double maxXdist = 0.;
    double maxYdist = 0.;
    for(auto connectedVertex : vertToVerts.at(v)) {
      Vertex *vj = mesh->getVertex(connectedVertex);
      maxXdist = fmax(maxXdist, 2.0 * fabs(vi->x() - vj->x()));
      maxYdist = fmax(maxYdist, 2.0 * fabs(vi->y() - vj->y()));
    }
    _scaling[v] = {maxXdist, maxYdist};
    // _scaling[v] = {fmax(maxXdist, maxYdist), fmax(maxXdist, maxYdist)};
  }

  // Some debug plots
  if(false) {
    system("rm vertexPatch*.pos");
    // Draw connected interior vertices
    for(auto v : vertToVerts) {
      // if(_isVertexBoundary[v.first]) {
      std::string posfile = "vertexPatch_P" + std::to_string(_elementOrder) + "_" + std::to_string(v.first) + ".pos";
      FILE *file = fopen(posfile.data(), "w");
      fprintf(file, "View\"vertexPatch%d\"{\n", v.first);
      Vertex *v1 = mesh->getVertex(v.first);
      fprintf(file, "SP(%g,%g,%g){%g};\n", v1->x(), v1->y(), 0., 2.);
      for(auto vertex : v.second) {
        Vertex *vv = mesh->getVertex(vertex);
        fprintf(file, "SP(%g,%g,%g){%g};\n", vv->x(), vv->y(), 0., 1.);
      }
      fprintf(file, "};\n");
      fclose(file);
      // }
    }
  }

  if(false) {
    system("rm interiorVertices*.pos");
    // Draw connected interior vertices
    for(auto &v : vertToInteriorVerts) {
      std::string posfile = "interiorVertices" + std::to_string(v.first) + ".pos";
      FILE *file = fopen(posfile.data(), "w");
      fprintf(file, "View\"interiorVertices\"{\n");
      Vertex *v1 = mesh->getVertex(v.first);
      fprintf(file, "SP(%g,%g,%g){%g};\n", v1->x(), v1->y(), 0., 2.);
      for(auto vertex : v.second) {
        Vertex *vv = mesh->getVertex(vertex);
        fprintf(file, "SP(%g,%g,%g){%g};\n", vv->x(), vv->y(), 0., 1.);
      }
      fprintf(file, "};\n");
      fclose(file);
    }
  }

  if(false) {
    system("rm elementPatch*.pos");
    // Draw patches of elements
    for(auto v : vertToElems) {
      std::string posfile = "elementPatch" + std::to_string(v.first) + ".pos";
      FILE *file = fopen(posfile.data(), "w");
      fprintf(file, "View\"patch\"{\n");
      for(auto elem : v.second) {
        mesh->getCoord(cnc->getTag(), elem, geoCoord);
        cnc->writeElementToPOS(file, geoCoord, 1.);
      }
      fprintf(file, "};\n");
      fclose(file);
    }
  }

   if(false) {
    system("rm scalingBox*.pos");
    for(auto pair : _scaling) {
      std::string posfile = "scalingBox_" + std::to_string(pair.first) + ".pos";
      FILE *file = fopen(posfile.data(), "w");
      fprintf(file, "View\"scalingBox_%d\"{\n", pair.first);
      Vertex *v1 = mesh->getVertex(pair.first);
      double lx = pair.second.first / 2.;
      double ly = pair.second.second / 2.;
      fprintf(file, "SL(%g,%g, 0.,%g,%g,0.){%g,%g};\n", v1->x() - lx, v1->y() - ly, v1->x() + lx, v1->y() - ly, 1., 1.);
      fprintf(file, "SL(%g,%g, 0.,%g,%g,0.){%g,%g};\n", v1->x() + lx, v1->y() - ly, v1->x() + lx, v1->y() + ly, 1., 1.);
      fprintf(file, "SL(%g,%g, 0.,%g,%g,0.){%g,%g};\n", v1->x() + lx, v1->y() + ly, v1->x() - lx, v1->y() + ly, 1., 1.);
      fprintf(file, "SL(%g,%g, 0.,%g,%g,0.){%g,%g};\n", v1->x() - lx, v1->y() + ly, v1->x() - lx, v1->y() - ly, 1., 1.);
      fprintf(file, "};\n");
      fclose(file);
    }
  }
}

//
// Increase the patch of elements and vertices around "vertex" by one layer
//
void feNewPatch::increasePatchSize(feMesh *mesh, const feCncGeo *cnc, const int vertex)
{
  if(isOnEdge.at(vertex)) {
    feInfo("Could not increase patch of edge vertex");
    exit(-1);
  }

  // Debug plot
  if(false) {
    // Draw connected interior vertices
    std::string posfile = "vertexPatch_beforeIncrease.pos";
    FILE *file = fopen(posfile.data(), "w");
    fprintf(file, "View\"vertexPatch%d\"{\n", vertex);
    Vertex *v1 = mesh->getVertex(vertex);
    fprintf(file, "SP(%g,%g,%g){%g};\n", v1->x(), v1->y(), 0., 2.);
    for(auto connectedVertex : vertToVerts.at(vertex)) {
      Vertex *vv = mesh->getVertex(connectedVertex);
      fprintf(file, "SP(%g,%g,%g){%g};\n", vv->x(), vv->y(), 0., 1.);
    }
    fprintf(file, "};\n");
    fclose(file);
  }

  if(_useOriginalZhangNagaPatchDefinition) {
    // For the original patch based on longest edge ball
    // Double the hMax and try again
    _hMax.at(vertex) *= 2.;
    double hMax = _hMax.at(vertex);

    // Loop and check for new vertices
    Vertex *vi = mesh->getVertex(vertex);

    for(auto vertex_j : _vertices) {
      Vertex *vj = mesh->getVertex(vertex_j);
      double dist =
        sqrt((vi->x() - vj->x()) * (vi->x() - vj->x()) + (vi->y() - vj->y()) * (vi->y() - vj->y()));
      if(dist <= 1.01 * hMax) {
        vertToVerts.at(vertex).insert(vertex_j);
      }
    }

  } else {
    // Add a layer of elements to the set in vertToElems
    const std::vector<int> &connecNodes = cnc->getVerticesConnectivity();
    std::set<int> toAdd;
    std::set<int> &connectedElements = vertToElems.at(vertex);

    for(auto currentElem : connectedElements) {
      for(int j = 0; j < _nNodePerElm; ++j) {
        for(auto elemToAdd : vertToElems_singleLayer[connecNodes[_nNodePerElm * currentElem + j]]) {
          toAdd.insert(elemToAdd);
        }
      }
    }

    for(auto e : toAdd) {
      connectedElements.insert(e);
    }

    // Add the vertices of the connected elements set to vertToVerts
    for(auto elem : connectedElements) {
      int nVerticesPerElem = _reconstructAtHighOrderNodes ? _nNodePerElm : 3;
      for(int j = 0; j < nVerticesPerElem; ++j) {
        vertToVerts[vertex].insert(connecNodes[_nNodePerElm * elem + j]);
      }
    }
  }

  // Debug plot
  if(false) {
    // Draw connected interior vertices
    std::string posfile = "vertexPatch_afterIncrease_" + std::to_string(vertex) + ".pos";
    FILE *file = fopen(posfile.data(), "w");
    fprintf(file, "View\"vertexPatch_afterIncrease%d\"{\n", vertex);
    Vertex *v1 = mesh->getVertex(vertex);
    fprintf(file, "SP(%g,%g,%g){%g};\n", v1->x(), v1->y(), 0., 2.);
    for(auto connectedVertex : vertToVerts.at(vertex)) {
      Vertex *vv = mesh->getVertex(connectedVertex);
      fprintf(file, "SP(%g,%g,%g){%g};\n", vv->x(), vv->y(), 0., 1.);
    }
    fprintf(file, "};\n");
    fclose(file);
  }
}

void feNewRecovery::setDimensions()
{
  // The dimension of the polynomial bases for recoveries is the one of degree k+1 :
  if(_dim == 0) {
    feErrorMsg(FE_STATUS_ERROR, "No recovery method for connectivity of dimension 0.");
    exit(-1);
  } else if(_dim == 1) {
    _dimRecovery = _degRec + 1;
    _dimDerivation = _degSol + 1;
  } else if(_dim == 2) {
    if(_cnc->getGeometry() == geometryType::TRI) {
      _dimRecovery = (_degRec + 1) * (_degRec + 2) / 2;
      _dimDerivation = (_degSol + 1) * (_degSol + 2) / 2;
    } else if(_cnc->getGeometry() == geometryType::QUAD) {
      _dimRecovery = (_degRec + 1) * (_degRec + 1);
      _dimDerivation = (_degSol + 1) * (_degSol + 1);
    } else {
      feErrorMsg(FE_STATUS_ERROR, "No polynomial recovery defined for geometry \"%s\".",
                 toString(_cnc->getGeometry()).data());
    }
  } else if(_dim == 3) {
    feErrorMsg(FE_STATUS_ERROR, "No polynomial recovery defined for geometry \"%s\".",
               toString(_cnc->getGeometry()).data());
  }
}

void feNewRecovery::setPolynomialExponents()
{
  // The exponents of the monomials :
  _expX_recov.resize(_dimRecovery, 0);
  _expY_recov.resize(_dimRecovery, 0);

  _expX_deriv.resize(_dimDerivation, 0);
  _expY_deriv.resize(_dimDerivation, 0);

  int ind = 0, n = _degRec;

  if(_dim == 1) {
    // The 1D basis is 1 x x^2 x^3 ...
    for(int i = 0; i < _dimRecovery; ++i) {
      _expX_recov[i] = i;
      if(i < _dimDerivation) _expX_deriv[i] = i;
    }
  } else if(_dim == 2) {
    ind = 0;
    for(int j = 0; j <= n; ++j) {
      for(int i = 0; i <= n - j; ++i) {
        _expX_recov[ind] = i;
        _expY_recov[ind] = j;
        ++ind;
      }
    }
    ind = 0;
    for(int j = 0; j < n; ++j) {
      for(int i = 0; i < n - j; ++i) {
        _expX_deriv[ind] = i;
        _expY_deriv[ind] = j;
        ++ind;
      }
    }

    // switch(_degRec) {
    //   case 1:
    //     _expX_recov = {0,1,0};
    //     _expY_recov = {0,0,1};
    //     _expX_deriv = {1};
    //     _expY_deriv = {1};
    //     break;
    //   case 2:
    //     _expX_recov = {0,1,0,2,1,0};
    //     _expY_recov = {0,0,1,0,1,2};
    //     _expX_deriv = {0,1,0};
    //     _expY_deriv = {0,0,1};
    //     break;
    //   case 3:
    //     _expX_recov = {0,1,0,2,1,0,3,2,1,0};
    //     _expY_recov = {0,0,1,0,1,2,0,1,2,3};
    //     _expX_deriv = {0,1,0,2,1,0};
    //     _expY_deriv = {0,0,1,0,1,2};
    //     break;
    //   case 4:
    //     _expX_recov = {0,1,0,2,1,0,3,2,1,0,4,3,2,1,0};
    //     _expY_recov = {0,0,1,0,1,2,0,1,2,3,0,1,2,3,4};
    //     _expX_deriv = {0,1,0,2,1,0,3,2,1,0};
    //     _expY_deriv = {0,0,1,0,1,2,0,1,2,3};
    //     break;
    //   default:
    //     feErrorMsg(FE_STATUS_ERROR, "No polynomial recovery defined for degree = %d.", _degRec);
    //     exit(-1);
    // }

  } else {
    feErrorMsg(FE_STATUS_ERROR, "No polynomial recovery defined for dim = %d.", _dim);
    exit(-1);
  }
}

void feNewRecovery::computeVertexMassMatrices2D()
{
  int nQuad = _geoSpace->getNumQuadPoints();
  const std::vector<double> &w = _geoSpace->getQuadratureWeights();
  const std::vector<double> &J = _cnc->getJacobians();

  std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  // Mass matrices centered at vertices
  const std::vector<int> vertices = _patch->getVertices();
#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
    std::vector<double> x(3, 0.0);
    std::vector<double> xLoc(3, 0.0);
    std::vector<double> monomials(_dimRecovery, 0.);

#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t iVertex = 0; iVertex < vertices.size(); ++iVertex) {
      int v = vertices[iVertex];
      double xv = _mesh->getVertex(v)->x();
      double yv = _mesh->getVertex(v)->y();

      // Get patch of elements
      const std::set<int> &elemPatch = _patch->getPatch(v);

      SquareMatrix massMatrix(_dimRecovery);

      for(auto elem : elemPatch) {
        _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

        // Loop over quad points and increment mass matrix
        for(int k = 0; k < nQuad; ++k) {
          double jac = J[nQuad * elem + k];

          // Get coordinates in the frame centered at vertex v
          _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
          xLoc[0] = x[0] - xv;
          xLoc[1] = x[1] - yv;
#if defined(APPLY_SCALING)
          xLoc[0] /= scaling.at(v).first;
          xLoc[1] /= scaling.at(v).second;
#endif

          for(int i = 0; i < _dimRecovery; ++i) {
            monomials[i] = pow(xLoc[0], _expX_recov[i]) * pow(xLoc[1], _expY_recov[i]);
          }

          for(int i = 0; i < _dimRecovery; ++i) {
            for(int j = 0; j < _dimRecovery; ++j) {
              massMatrix(i, j) += jac * w[k] * monomials[i] * monomials[j];
            }
          }
        }
      }

// We should try initializing the map beforehand and remove the critical,
// there seems to be data race when it's not initialized
#if defined(HAVE_OMP)
#pragma omp critical
#endif
      {
        _inverseMassSquareMatrices[v] = massMatrix.inverse();
      }
    }
  }
}

void feNewRecovery::computeRHSAndSolve(int numRecoveries, int nRecoveredFields,
                                       int nRecoveredDerivatives, int iDerivative)
{
  int nQuad = _geoSpace->getNumQuadPoints();
  const std::vector<double> &w = _geoSpace->getQuadratureWeights();
  const std::vector<double> &J = _cnc->getJacobians();

  std::vector<int> &vertices = _patch->getVertices();
  std::map<int, bool> &isVertexBoundary = _patch->getBoundaryVertices();
  std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  int _nVertPerElm = _nNodePerElm;

  std::vector<double> &solVec = _sol->getSolutionReference();

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    std::vector<feInt> adr(_intSpace->getNumFunctions());
    std::vector<double> solution(_intSpace->getNumFunctions());
    std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
    std::vector<double> x(3, 0.0);
    std::vector<double> xLoc(3, 0.0);
    std::vector<double> monomials(_dimRecovery, 0.);
    std::vector<double> recoveryVector(_dimRecovery, 0.);
    std::vector<double> recoveryVectorOnBoundary(1, 0.);
    std::vector<double> u(numRecoveries, 0.);

#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(auto v : vertices) {
#if defined(TREAT_BOUNDARIES)
      // Only reconstruct on non-boundary vertices
      if(!isVertexBoundary[v]) {
#endif
        double xv = _mesh->getVertex(v)->x();
        double yv = _mesh->getVertex(v)->y();

        std::vector<Vector> RHS;
        RHS.clear();
        for(int i = 0; i < numRecoveries; ++i) RHS.push_back(Vector(_dimRecovery));

        // Get patch of elements
        std::set<int> &elemPatch = _patch->getPatch(v);

        for(auto elem : elemPatch) {
          _intSpace->initializeAddressingVector(elem, adr);
          for(size_t i = 0; i < adr.size(); ++i) {
            solution[i] = solVec[adr[i]];
          }
          _mesh->getCoord(_cnc, elem, geoCoord);

          // Loop over quad points and increment right hand side
          for(int k = 0; k < nQuad; ++k) {
            double jac = J[nQuad * elem + k];

            _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
            xLoc[0] = x[0] - xv;
            xLoc[1] = x[1] - yv;
#if defined(APPLY_SCALING)
            xLoc[0] /= scaling.at(v).first;
            xLoc[1] /= scaling.at(v).second;
#endif

            for(int i = 0; i < _dimRecovery; ++i) {
              monomials[i] = pow(xLoc[0], _expX_recov[i]) * pow(xLoc[1], _expY_recov[i]);
            }

            if(iDerivative == 0) {
              // Reconstruct the solution : Simply interpolate the solution at quad nodes
              for(int iComp = 0; iComp < numRecoveries; ++iComp) {
                u[iComp] = _intSpace->interpolateFieldAtQuadNode(solution, k);
              }
            } else {
              // Reconstruct the derivatives
              // The contributions to the derivative from the vertices must be
              // evaluated and averaged at quad nodes to avoid a trivial solution
              for(int iComp = 0; iComp < numRecoveries; ++iComp) {
                // Average the evaluations at edge vertices, then interpolate using
                // Lagrange functions (Zhang & Naga)
                u[iComp] = evaluateRecoveryAtQuadNode(PPR::DERIVATIVE,
                                                      nRecoveredDerivatives + iComp, elem, k);
              }
            }

            for(int iComp = 0; iComp < numRecoveries; ++iComp) {
              for(int i = 0; i < _dimRecovery; ++i) {
                RHS[iComp](i) += jac * w[k] * u[iComp] * monomials[i];
              }
            }
          }
        }

        for(int iComp = 0; iComp < numRecoveries; ++iComp) {
          Vector sol = _inverseMassSquareMatrices[v] * RHS[iComp];

          for(int i = 0; i < _dimRecovery; ++i) {
            recoveryVector[i] = sol(i);
          }

#if defined(HAVE_OMP)
#pragma omp critical
#endif
          recoveryCoeff[v][nRecoveredFields + iComp] = recoveryVector;
        }
#if defined(TREAT_BOUNDARIES)
      } // if not a boundary vertex
#endif
    }

#if defined(TREAT_BOUNDARIES)
// Wait for all recoveries to be computed at interior vertices
#if defined(HAVE_OMP)
#pragma omp barrier
#endif

// Now treat the boundary vertices
#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(auto v : vertices) {
      // Only reconstruct on non-boundary vertices
      if(isVertexBoundary[v]) {
        // Get patch of interior vertices
        std::set<int> &connectedInteriorVertices = _patch->getPatchOfInteriorVertices(v);

        for(int iComp = 0; iComp < numRecoveries; ++iComp) {
          recoveryVectorOnBoundary[0] = evaluatePolynomialAtBoundaryVertex(
            PPR::RECOVERY, nRecoveredFields + iComp, v, connectedInteriorVertices);

#if defined(HAVE_OMP)
#pragma omp critical
#endif
          recoveryCoeff[v][nRecoveredFields + iComp] = recoveryVectorOnBoundary;
        }
      } // if boundary vertex
    }
#endif
  } // pragma omp parallel
}

void feNewRecovery::computeVertexLeastSquareMatrices2D()
{
  const std::vector<int> &vertices = _patch->getVertices();
  const std::map<int, std::set<int> > &vertToVerts = _patch->getVerticesPatches();
  const std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  int maxPatchIncreases = 5;

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double xLoc[2];

#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t iVertex = 0; iVertex < vertices.size(); ++iVertex) {
      int v = vertices[iVertex];
      double xv = _mesh->getVertex(v)->x();
      double yv = _mesh->getVertex(v)->y();

      const std::set<int> &verticesPatch = vertToVerts.at(v);

      bool rankIsOk = false;
      bool rankWasBad = false;
      int numPatchIncrease = 0;

      do {
        int numConnectedVertices = verticesPatch.size();
        if(numConnectedVertices < _dimRecovery) {
          feInfo(
            "Not enough connected vertices (%d) - This should not happen because the number of "
            "connected vertices should be at least the dimension of the recovery basis (rank is "
            "checked after)",
            numConnectedVertices);
          exit(-1);
        }

        Matrix A(numConnectedVertices, _dimRecovery);

        int iV = 0;
        for(auto patchVertex : verticesPatch) {
          // Get coordinates in the frame centered at vertex v
          double xvj = _mesh->getVertex(patchVertex)->x();
          double yvj = _mesh->getVertex(patchVertex)->y();
          xLoc[0] = xvj - xv;
          xLoc[1] = yvj - yv;
#if defined(APPLY_SCALING)
          xLoc[0] /= scaling.at(v).first;
          xLoc[1] /= scaling.at(v).second;
#endif
          // if(fabs(xLoc[0]) < XLOC_TOL) xLoc[0] = 0.;
          // if(fabs(xLoc[1]) < XLOC_TOL) xLoc[1] = 0.;

          for(int j = 0; j < _dimRecovery; ++j) {
            A(iV, j) = myPow(xLoc[0], _expX_recov[j]) * myPow(xLoc[1], _expY_recov[j]);
          }
          iV++;
        }

        int rankA = A.rank();
        SquareMatrix foo = A.multiplyByTranspose();
        int rankAtA = foo.rank();

        if(fmin(rankA, rankAtA) >= _dimRecovery) {
          rankIsOk = true;

// Get (A^T*A)^-1 * A^T
#if defined(HAVE_OMP)
#pragma omp critical
#endif
          _leastSquaresMatrices[v] = A.getLeastSquaresMatrix();

        } else {
          if(numPatchIncrease++ > maxPatchIncreases) {
            feInfo(
              "Cannot compute recovery at vertex %d because least square matrix rank is = %d < %d",
              v, fmin(rankA, rankAtA), _dimRecovery);
            feInfo("Could not satisfy rank condition even after %d increases in the patch size!",
                   maxPatchIncreases);
            exit(-1);
          }

#if defined(HAVE_OMP)
// Critical because modifying maps in the patch structure
#pragma omp critical
#endif
          _patch->increasePatchSize(_mesh, _cnc, v);
        }
      } while(!rankIsOk);
    }
  }
}

void feNewRecovery::recomputeLeastSquareMatrix(const int vertex)
{
  const std::map<int, std::set<int> > &vertToVerts = _patch->getVerticesPatches();
  const std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  double xLoc[2];
  double xv = _mesh->getVertex(vertex)->x();
  double yv = _mesh->getVertex(vertex)->y();

  const std::set<int> &verticesPatch = vertToVerts.at(vertex);

  int numConnectedVertices = verticesPatch.size();
  if(numConnectedVertices < _dimRecovery) {
    feInfo("Not enough connected vertices (%d)", numConnectedVertices);
    exit(-1);
  }

  Matrix A(numConnectedVertices, _dimRecovery);

  int iV = 0;
  for(auto patchVertex : verticesPatch) {
    // Get coordinates in the frame centered at vertex v
    double xvj = _mesh->getVertex(patchVertex)->x();
    double yvj = _mesh->getVertex(patchVertex)->y();
    xLoc[0] = xvj - xv;
    xLoc[1] = yvj - yv;
#if defined(APPLY_SCALING)
    xLoc[0] /= scaling.at(vertex).first;
    xLoc[1] /= scaling.at(vertex).second;
#endif
    // if(fabs(xLoc[0]) < XLOC_TOL) xLoc[0] = 0.;
    // if(fabs(xLoc[1]) < XLOC_TOL) xLoc[1] = 0.;

    for(int j = 0; j < _dimRecovery; ++j) {
      A(iV, j) = myPow(xLoc[0], _expX_recov[j]) * myPow(xLoc[1], _expY_recov[j]);
    }
    iV++;
  }

  int rank = A.rank();
  if(rank < _dimRecovery) {
    feInfo("Cannot compute recovery at vertex %d because least square matrix rank is = %d < %d",
           vertex, rank, _dimRecovery);
    feInfo("Adjust the patch around this vertex.");
    exit(-1);
  }

  // Get (A^T*A)^-1 * A^T
  _leastSquaresMatrices[vertex] = A.getLeastSquaresMatrix();
}

void feNewRecovery::computeRHSAndSolve_noIntegral(int numRecoveries, int nRecoveredFields,
                                                  int nRecoveredDerivatives, int iDerivative)
{
  std::vector<int> &vertices = _patch->getVertices();
  std::map<int, std::set<int> > &vertToVerts = _patch->getVerticesPatches();
  std::map<int, bool> &isVertexBoundary = _patch->getBoundaryVertices();
  std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();
  std::vector<double> &solVec = _sol->getSolutionReference();

  std::map<int, Vector> allSOL;
  std::map<int, Vector> allRHS;

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    int elm;
    double r[3];
    // std::vector<double> xLoc(3, 0.0);
    std::vector<double> pos(3, 0.0);
    std::vector<feInt> adr(_intSpace->getNumFunctions());
    std::vector<double> sol(_intSpace->getNumFunctions());
    std::vector<double> recoveryVector(_dimRecovery, 0.);
    std::vector<double> recoveryVectorOnBoundary(1, 0.);
    std::vector<double> u(numRecoveries, 0.);

#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t iVertex = 0; iVertex < vertices.size(); ++iVertex) {
      int v = vertices[iVertex];
      std::set<int> &verticesPatch = vertToVerts[v];
      int numConnectedVertices = verticesPatch.size();

      double xv = _mesh->getVertex(v)->x();
      double yv = _mesh->getVertex(v)->y();

      std::vector<Vector> RHS;
      RHS.clear();
      for(int i = 0; i < numRecoveries; ++i) RHS.push_back(Vector(numConnectedVertices));

      int iV = 0;
      for(auto patchVertex : verticesPatch) {
        double xvj = _mesh->getVertex(patchVertex)->x();
        double yvj = _mesh->getVertex(patchVertex)->y();
//         xLoc[0] = xvj - xv;
//         xLoc[1] = yvj - yv;
// #if defined(APPLY_SCALING)
//         xLoc[0] /= scaling.at(v).first;
//         xLoc[1] /= scaling.at(v).second;
// #endif
        // if(fabs(xLoc[0]) < XLOC_TOL) xLoc[0] = 0.;
        // if(fabs(xLoc[1]) < XLOC_TOL) xLoc[1] = 0.;

        if(_numbering == nullptr) {
          feInfo("Could not compute pointwise least squares for solution recovery because "
                 "numbering is NULL");
          exit(-1);
        }

        if(iDerivative == 0) {
          // Reconstruct the solution
          int iDOF = _numbering->getNumbering(_intSpace->getFieldID())
                       ->getDOFNumberAtVertex(patchVertex, _componentToRecover);
          for(int iComp = 0; iComp < numRecoveries; ++iComp) {
            if(iDOF >= 0) {
              // There is a degree of freedom at this mesh vertex
              // Get solution from feSolution at DOF
              RHS[iComp](iV) = _sol->getSolAtDOF(iDOF);
            } else {
              // No DOF associated to the mesh vertex, e.g. P2 mesh vertex
              // and P3 interpolant. Interpolate solution at vertex.
              pos[0] = xvj;
              pos[1] = yvj;
              // TODO: Search only in the element patch or better directly provide the reference
              // coord if there is a way to obtain them easily
              bool wasFound = _mesh->locateVertex(pos.data(), elm, r, 1e-6);
              if(!wasFound) {
                feInfo(
                  "Could not locate point %f - %f - %f in the mesh when reconstructing solution!",
                  pos[0], pos[1], pos[2]);
                exit(-1);
              }
              _intSpace->initializeAddressingVector(elm, adr);

              for(size_t i = 0; i < adr.size(); ++i) {
                sol[i] = solVec[adr[i]];
              }

              RHS[iComp](iV) =
                _intSpace->interpolateVectorFieldComponent(sol, _componentToRecover, r);
            }
          }
        } else {
          // Reconstruct the derivatives
          for(int iComp = 0; iComp < numRecoveries; ++iComp) {
            RHS[iComp](iV) = derivativeCoeff[patchVertex][nRecoveredDerivatives + iComp][0];
          }
        }
        ++iV;
      }

      for(int iComp = 0; iComp < numRecoveries; ++iComp) {

        Vector sol = _leastSquaresMatrices[v] * RHS[iComp];

        #if defined(APPLY_SCALING)
        // Divide solution vector by scaling coefficient with the right exponent
        for(int iii = 0; iii < _dimRecovery; ++iii) {
          double h = pow(scaling.at(v).first, _expX_recov[iii]) * pow(scaling.at(v).second, _expY_recov[iii]);
         sol(iii) = sol(iii) / h;
        }
        #endif

        for(int i = 0; i < _dimRecovery; ++i) {
          recoveryVector[i] = sol(i);
        }

#if defined(HAVE_OMP)
#pragma omp critical
#endif
        {
          recoveryCoeff[v][nRecoveredFields + iComp] = recoveryVector;
        }
      }
    }
  } // pragma omp parallel
}

// Solve the least squares system at a given vertex and for each component of the input metric
// Used to recover the derivatives of M^(-1/2) for curved metrics
void feNewRecovery::computeRHSAndSolve_noIntegral_inputMetric(
  const std::map<int, int> &patchVertices2nodeTags,
  const std::map<int, MetricTensor> &inputMetrics,
  std::map<int, MetricTensor> &dMdx,
  std::map<int, MetricTensor> &dMdy)
{
  std::vector<int> &vertices = _patch->getVertices();
  std::map<int, std::set<int> > &vertToVerts = _patch->getVerticesPatches();
  std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  // Number of independent components of the symmetric 2D metric tensors
  int numRecoveries = 3; 

  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    std::vector<double> sol(_intSpace->getNumFunctions());
    std::vector<double> recoveryVector(_dimRecovery, 0.);

    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t iVertex = 0; iVertex < vertices.size(); ++iVertex) {
      int v = vertices[iVertex];
      std::set<int> &verticesPatch = vertToVerts[v];
      int numConnectedVertices = verticesPatch.size();

      // Create and fill the three RHS vectors
      std::vector<Vector> RHS;
      RHS.clear();
      for(int i = 0; i < numRecoveries; ++i) {
        RHS.push_back(Vector(numConnectedVertices));
      }

      int iV = 0;
      for(auto patchVertex : verticesPatch) {
        int ind = patchVertices2nodeTags.at(patchVertex);
        double m11 = inputMetrics.at( ind )(0,0);
        double m12 = inputMetrics.at( ind )(0,1);
        double m22 = inputMetrics.at( ind )(1,1);
        RHS[0](iV) = m11;
        RHS[1](iV) = m12;
        RHS[2](iV) = m22;
        ++iV;
      }

      // Compute least-squares solution
      for(int iComp = 0; iComp < numRecoveries; ++iComp) {

        Vector sol = _leastSquaresMatrices[v] * RHS[iComp];

        for(int i = 0; i < _dimRecovery; ++i) {
          #if defined(APPLY_SCALING)
          // Divide solution vector by scaling coefficient with the right exponent
          double h = pow(scaling.at(v).first, _expX_recov[i]) * pow(scaling.at(v).second, _expY_recov[i]);
          sol(i) = sol(i) / h;
          #endif
          recoveryVector[i] = sol(i);
        }

        // Compute the derivative at (0,0)
        double dmdx, dmdy;

        // Compute only the coefficient at 0
        for(int i = 0; i < _dimRecovery; ++i) {
          if(_expX_recov[i] == 1 && _expY_recov[i] == 0) {
            dmdx = recoveryVector[i];
            break;
          }
        }
        for(int i = 0; i < _dimRecovery; ++i) {
          if(_expX_recov[i] == 0 && _expY_recov[i] == 1) {
            dmdy = recoveryVector[i];
            break;
          }
        }

        int indv = patchVertices2nodeTags.at(v);

        if(iComp == 0) {
          dMdx.at( indv )(0,0) = dmdx;
          dMdy.at( indv )(0,0) = dmdy;
        } else if(iComp == 1) {
          dMdx.at( indv )(0,1) = dmdx;
          dMdy.at( indv )(0,1) = dmdy;
          dMdx.at( indv )(1,0) = dmdx;
          dMdy.at( indv )(1,0) = dmdy;
        } else {
          dMdx.at( indv )(1,1) = dmdx;
          dMdy.at( indv )(1,1) = dmdy;
        }
        
      }
    }
  } // pragma omp parallel
}

// Recompute the solution recovery at a single vertex (used when curving to minimize interpolation
// error)
void feNewRecovery::recomputeRHSAndSolve(const int vertex)
{
  std::vector<double> recoveryVector(_dimRecovery, 0.);
  std::map<int, std::set<int> > &vertToVerts = _patch->getVerticesPatches();
  std::set<int> &verticesPatch = vertToVerts[vertex];
  int numConnectedVertices = verticesPatch.size();

  Vector RHS(numConnectedVertices);

  int iV = 0;
  for(auto patchVertex : verticesPatch) {
    // Reconstruct the solution : get solution from feSolution at DOF
    int iDOF = _numbering->getNumbering(_intSpace->getFieldID())
                 ->getDOFNumberAtVertex(patchVertex, _componentToRecover);
    RHS(iV) = _sol->getSolAtDOF(iDOF);
    ++iV;
  }

  Vector sol = _leastSquaresMatrices[vertex] * RHS;

  for(int i = 0; i < _dimRecovery; ++i) {
    recoveryVector[i] = sol(i);
  }

  recoveryCoeff[vertex][0] = recoveryVector;
  recoveryIndependantTerm[_numTotalRecoveries * vertex + 0] = recoveryVector[0];
}

void feNewRecovery::computeDerivative(int nRecoveredFields, int indexRecovery, int iDerivative,
                                      std::ostream &output)
{
  std::vector<int> &vertices = _patch->getVertices();
  std::map<int, bool> &isVertexBoundary = _patch->getBoundaryVertices();
  const std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();

  // If a coefficient in the derivative is less than this relative threshold,
  // it is ignored
  double tol = 0.;

  if(_dim == 1) {
    for(auto v : vertices) {
      std::vector<double> &u = recoveryCoeff[v][nRecoveredFields + indexRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      for(int i = 0; i < _dimDerivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeff[v][indexRecovery] = dudx;
    }

  } else if(_dim == 2) {
    for(auto v : vertices) {
#if defined(TREAT_BOUNDARIES)
      if(!isVertexBoundary[v]) {
#endif
        std::vector<double> &u = recoveryCoeff[v][nRecoveredFields + indexRecovery];
        int indX = 0, indY = 0;
        std::vector<double> dudx(_dimDerivation, 0.);
        std::vector<double> dudy(_dimDerivation, 0.);

        double sumCoeffX = 0.0;
        double sumCoeffY = 0.0;

        for(int i = 0; i < _dimRecovery; ++i) {
          if(_expX_recov[i] != 0) {
            sumCoeffX += ((double)_expX_recov[i]) * fabs(u[i]);
          }
          if(_expY_recov[i] != 0) {
            sumCoeffY += ((double)_expY_recov[i]) * fabs(u[i]);
          }
        }

        for(int i = 0; i < _dimRecovery; ++i) {
          if(_expX_recov[i] != 0) {
            if(fabs(u[i] / sumCoeffX) > tol) {
              dudx[indX] = _expX_recov[i] * u[i];
            }
            indX++;
          }
          if(_expY_recov[i] != 0) {
            if(fabs(u[i] / sumCoeffY) > tol) {
              dudy[indY] = _expY_recov[i] * u[i];
            }
            indY++;
          }
        }

        derivativeCoeff[v][2 * nRecoveredFields + 2 * indexRecovery + 0] = dudx;
        derivativeCoeff[v][2 * nRecoveredFields + 2 * indexRecovery + 1] = dudy;
#if defined(TREAT_BOUNDARIES)
      }
#endif
    }

#if defined(TREAT_BOUNDARIES)
    // Evaluate derivative at boundary vertices
    for(auto v : vertices) {
      if(isVertexBoundary[v]) {
        std::set<int> &connectedInteriorVertices = _patch->getPatchOfInteriorVertices(v);
        std::vector<double> dudx(1), dudy(1);
        int index_dudx = 2 * nRecoveredFields + 2 * indexRecovery + 0;
        int index_dudy = 2 * nRecoveredFields + 2 * indexRecovery + 1;
        dudx[0] = evaluatePolynomialAtBoundaryVertex(PPR::DERIVATIVE, index_dudx, v,
                                                     connectedInteriorVertices);
        dudy[0] = evaluatePolynomialAtBoundaryVertex(PPR::DERIVATIVE, index_dudy, v,
                                                     connectedInteriorVertices);
        derivativeCoeff[v][index_dudx] = dudx;
        derivativeCoeff[v][index_dudy] = dudy;
      }
    }
#endif
  }

#if defined(EXPORT_FIELDS_TO_GMSH)
  // Print the solution at first pass
  if(iDerivative == 0) {
    output << "$NodeData\n";
    output << "1\n\"" << _intSpace->getFieldID() << "\"\n1\n3000\n3\n0\n1\n"
           << vertices.size() << "\n";
    for(auto v : vertices) {
      output << _mesh->getVertex(v)->getTag() << " " << std::setprecision(16)
             << recoveryCoeff.at(v)[indexRecovery][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }

  // Print solution at each pass
  // output << "$NodeData\n";
  // output << "1\n\"" << _intSpace->getFieldID() << "\"\n1\n3000\n3\n0\n1\n"
  //        << vertices.size() << "\n";
  // for(auto v : vertices) {
  //   output << _mesh->getVertex(v)->getTag() << " " << recoveryCoeff[v][nRecoveredFields +
  //   indexRecovery][0]
  //          << std::endl;
  // }
  // output << "$EndNodeData\n";

  // Also write the derivatives to the common file
  for(int i = 0; i < _dim; ++i) {
    output << "$NodeData\n";
    std::string fieldName = "d" + std::to_string(iDerivative + 1) + _intSpace->getFieldID() +
                            suffix[{iDerivative + 1, _dim * indexRecovery + i}];
    output << "1\n\"" << fieldName << "\"\n1\n3000\n3\n0\n1\n" << vertices.size() << "\n";
    for(auto v : vertices) {
      output << _mesh->getVertex(v)->getTag() << " " << std::setprecision(16)
             << derivativeCoeff.at(v)[_dim * nRecoveredFields + _dim * indexRecovery + i][0]
             << std::endl;
    }
    output << "$EndNodeData\n";
  }
#endif
}

thread_local std::vector<double> DUDX(10, 0.);
thread_local std::vector<double> DUDY(10, 0.);

void feNewRecovery::recomputeFirstDerivative(const int vertex)
{
  // If a coefficient in the derivative is less than this relative threshold,
  // it is ignored
  double tol = 0.;

  std::vector<double> &u = recoveryCoeff[vertex][0];
  int indX = 0, indY = 0;

  double sumCoeffX = 0.0;
  double sumCoeffY = 0.0;

  for(int i = 0; i < _dimRecovery; ++i) {
    if(_expX_recov[i] != 0) {
      sumCoeffX += ((double)_expX_recov[i]) * fabs(u[i]);
    }
    if(_expY_recov[i] != 0) {
      sumCoeffY += ((double)_expY_recov[i]) * fabs(u[i]);
    }
  }

  for(int i = 0; i < _dimRecovery; ++i) {
    if(_expX_recov[i] != 0) {
      if(fabs(u[i] / sumCoeffX) > tol) {
        DUDX[indX] = _expX_recov[i] * u[i];
      }
      indX++;
    }
    if(_expY_recov[i] != 0) {
      if(fabs(u[i] / sumCoeffY) > tol) {
        DUDY[indY] = _expY_recov[i] * u[i];
      }
      indY++;
    }
  }

  derivativeCoeff[vertex][0] = DUDX;
  derivativeCoeff[vertex][1] = DUDY;
  derivativeIndependantTerm[_numTotalDerivatives * vertex + 0] = DUDX[0];
  derivativeIndependantTerm[_numTotalDerivatives * vertex + 1] = DUDY[0];
}

void numYExponents(int pos, int nLoc, int &ny)
{
  if(nLoc == 1) {
    return;
  } else {
    if(pos >= nLoc / 2) {
      ny++;
      pos -= nLoc / 2;
    }
    numYExponents(pos, nLoc / 2, ny);
  }
}

void feNewRecovery::computeHomogeneousErrorPolynomials()
{
  if(_dim == 2) {
    std::vector<int> &vertices = _patch->getVertices();
    std::map<int, bool> &isVertexBoundary = _patch->getBoundaryVertices();

    std::vector<double> error(_degSol + 2, 0.);

    // Number of derivatives of the highest order
    int nDerivatives = pow(2, _degSol + 1);
    std::vector<int> target(nDerivatives);

    // High order derivatives contributes to the
    // entry in the error vector according to
    // its y exponent:
    // Example for cubic homogeneous polynomial:
    // error[0] = x^3 * y^0 coefficient
    // error[1] = x^2 * y^1 coefficient
    // error[2] = x^1 * y^2 coefficient
    // error[3] = x^0 * y^3 coefficient
    for(int i = 0; i < nDerivatives; ++i) {
      int ny = 0;
      numYExponents(i, nDerivatives, ny);
      target[i] = ny;
    }

    // Number of computed derivatives of lower order
    int offset = pow(2, _degSol + 2) - 2 - nDerivatives;

    for(auto v : vertices) {
      std::fill(error.begin(), error.end(), 0.);

#if defined(TREAT_BOUNDARIES)
      if(isVertexBoundary[v]) {
        std::set<int> &connectedInteriorVertices = _patch->getPatchOfInteriorVertices(v);
        // Interior polynomials averaged at boundary vertex
        for(int i = 0; i < nDerivatives; ++i) {
          error[target[i]] += evaluatePolynomialAtBoundaryVertex(PPR::DERIVATIVE, offset + i, v,
                                                                 connectedInteriorVertices);
        }

      } else {
#endif
        // Interior vertex's own polynomial
        for(int i = 0; i < nDerivatives; ++i) {
          error[target[i]] += derivativeCoeff[v][offset + i][0];
        }
#if defined(TREAT_BOUNDARIES)
      }
#endif

      // Set minimum error to threshold
      for(auto &val : error) {
        if(val >= 0) val = fmax(1e-10, val);
        if(val <= 0) val = fmin(-1e-10, val);
      }

      errorCoeff[v] = error;
    }
  }
}

double feNewRecovery::evaluatePolynomial(PPR recoveredField, const int index, const int vertex,
                                         const double *xLoc)
{
  double res = 0.;
  if(recoveredField == PPR::RECOVERY) {
    // Evaluate recovery
    for(int i = 0; i < _dimRecovery; ++i) {
      res += myPow(xLoc[0], _expX_recov[i]) * myPow(xLoc[1], _expY_recov[i]) *
             recoveryCoeff.at(vertex)[index][i];
    }
  } else {
    // Evaluate derivative
    for(int i = 0; i < _dimDerivation; ++i) {
      res += myPow(xLoc[0], _expX_deriv[i]) * myPow(xLoc[1], _expY_deriv[i]) *
             derivativeCoeff.at(vertex)[index][i];
    }
  }
  return res;
}

// Evaluate polynomial at boundary vertex as the average of
// the evaluations of polynomials of connected interior vertices
double
feNewRecovery::evaluatePolynomialAtBoundaryVertex(PPR recoveredField, const int index,
                                                  const int vertex,
                                                  const std::set<int> &connectedInteriorVertices)
{
  double xLoc[2];
  std::map<int, std::pair<double, double> > &scaling = _patch->getScaling();
  Vertex *vInt, *vBnd = _mesh->getVertex(vertex);
  double res = 0.;
  for(auto v : connectedInteriorVertices) {
    vInt = _mesh->getVertex(v);
    xLoc[0] = vBnd->x() - vInt->x();
    xLoc[1] = vBnd->y() - vInt->y();
    // if(recoveredField == PPR::RECOVERY)
    //   feInfo("Evaluating at vertex %d - boundary = %d - edge = %d - size = %d", v,
    //   _patch->isVertexBoundary(v), _patch->isVertexOnEdge(v), recoveryCoeff[v][index].size());
    // else
    //   feInfo("Evaluating at vertex %d - boundary = %d - edge = %d - size = %d", v,
    //   _patch->isVertexBoundary(v), _patch->isVertexOnEdge(v), derivativeCoeff[v][index].size());
    res += evaluatePolynomial(recoveredField, index, v, xLoc);
  }
  if(connectedInteriorVertices.empty()) {
    feInfo("No connected interior vertex...");
    exit(-1);
  }
  res /= (double)connectedInteriorVertices.size();
  return res;
}

double feNewRecovery::evaluateRecoveryAtVertex(PPR recoveredField, const int index,
                                               const int vertex)
{
  if(recoveredField == PPR::RECOVERY) {
    return recoveryIndependantTerm[_numTotalRecoveries * vertex + index];
  } else {
    return derivativeIndependantTerm[_numTotalDerivatives * vertex + index];
  }
}

void feNewRecovery::computeTransformedEdgeLength(std::vector<double> &geoCoord,
                                                 std::vector<double> &dxdr,
                                                 std::vector<double> &dxds, const double tEnd,
                                                 int whichEdge, double &firstHalf,
                                                 double &secondHalf)
{
  double w1, w2, t1, t2, l1 = 0., l2 = 0.;
  double R[3] = {0., 0., 0.}, r[3], s[3];
  double drdt[3] = {1., -1., 0.};
  double dsdt[3] = {0., 1., -1.};

  for(size_t i = 0; i < _nQuad1d; ++i) {
    // Map from [-1,1] to parameter space t in [0,tEnd] or in [tEnd,1]
    t1 = _xQuad1d[i] * tEnd / 2. + tEnd / 2.;
    w1 = _wQuad1d[i] * tEnd / 2.;
    t2 = _xQuad1d[i] * (1. - tEnd) / 2. + (1. + tEnd) / 2.;
    w2 = _wQuad1d[i] * (1. - tEnd) / 2.;

    r[0] = t1;
    r[1] = 1 - t1;
    r[2] = 0.;

    s[0] = 0.;
    s[1] = t1;
    s[2] = 1. - t1;

    R[0] = r[whichEdge];
    R[1] = s[whichEdge];
    R[2] = 0.;

    _geoSpace->interpolateVectorField_rDerivative(geoCoord, R, dxdr);
    _geoSpace->interpolateVectorField_sDerivative(geoCoord, R, dxds);

    double dxdt = dxdr[0] * drdt[whichEdge] + dxds[0] * dsdt[whichEdge];
    double dydt = dxdr[1] * drdt[whichEdge] + dxds[1] * dsdt[whichEdge];
    l1 += w1 * sqrt(dxdt * dxdt + dydt * dydt);

    r[0] = t2;
    r[1] = 1 - t2;
    r[2] = 0.;

    s[0] = 0.;
    s[1] = t2;
    s[2] = 1. - t2;

    R[0] = r[whichEdge];
    R[1] = s[whichEdge];
    R[2] = 0.;

    _geoSpace->interpolateVectorField_rDerivative(geoCoord, R, dxdr);
    _geoSpace->interpolateVectorField_sDerivative(geoCoord, R, dxds);

    dxdt = dxdr[0] * drdt[whichEdge] + dxds[0] * dsdt[whichEdge];
    dydt = dxdr[1] * drdt[whichEdge] + dxds[1] * dsdt[whichEdge];
    l2 += w2 * sqrt(dxdt * dxdt + dydt * dydt);
  }

  firstHalf = l1;
  secondHalf = l2;
}

double computeApproximateLength(feSpace *geoSpace, std::vector<double> &geoCoord,
                                                 int whichEdge,
                                                 double x0[2],
                                                 int N, double min_t, double max_t)
{
  double R[3] = {0., 0., 0.}, r[3], s[3];
  std::vector<double> pos(3, 0.);

  FILE *file = fopen("approxLength.pos", "w");
  fprintf(file, "View\"length\"{\n");

  double l = 0.; 
  double x[2] = {x0[0], x0[1]};
  double xprev[2] = {x[0], x[1]};
  double dt = (max_t - min_t) / N;
  double t = min_t;
  for(size_t i = 0; i < N+1; ++i) {

    fprintf(file, "SP(%1.14e,%1.14e,0.){%f};\n", x[0], x[1], 1.);
    
    r[0] = t;
    r[1] = 1 - t;
    r[2] = 0.;

    s[0] = 0.;
    s[1] = t;
    s[2] = 1. - t;

    R[0] = r[whichEdge];
    R[1] = s[whichEdge];
    R[2] = 0.;

    geoSpace->interpolateVectorField(geoCoord, R, pos);

    x[0] = pos[0];
    x[1] = pos[1];
    double dist = sqrt( (x[0] - xprev[0])*(x[0] - xprev[0]) + (x[1] - xprev[1])*(x[1] - xprev[1]) );
    l += dist;
    t += dt;
    xprev[0] = x[0];
    xprev[1] = x[1];
  }

  fprintf(file, "};\n");
  fclose(file);

  feInfo("ending with t = %f - t-dt = %f", t, t - dt);

  return l;
}

double feNewRecovery::computeTransformedEdgeLength2(std::vector<double> &geoCoord,
                                                 std::vector<double> &dxdr,
                                                 std::vector<double> &dxds,
                                                 const double t0,
                                                 const double tEnd,
                                                 int whichEdge)
{
  double w, t, l = 0.;
  double R[3] = {0., 0., 0.}, r[3], s[3];
  double drdt[3] = {1., -1., 0.};
  double dsdt[3] = {0., 1., -1.};

  for(size_t i = 0; i < _nQuad1d; ++i) {
    // Map from [-1,1] to parameter space t in [t0,tEnd]
    t = (tEnd - t0)/2. * _xQuad1d[i] + (tEnd + t0)/2.;
    w = _wQuad1d[i] * (tEnd - t0) / 2.;

    r[0] = t;
    r[1] = 1 - t;
    r[2] = 0.;

    s[0] = 0.;
    s[1] = t;
    s[2] = 1. - t;

    R[0] = r[whichEdge];
    R[1] = s[whichEdge];
    R[2] = 0.;

    _geoSpace->interpolateVectorField_rDerivative(geoCoord, R, dxdr);
    _geoSpace->interpolateVectorField_sDerivative(geoCoord, R, dxds);

    double dxdt = dxdr[0] * drdt[whichEdge] + dxds[0] * dsdt[whichEdge];
    double dydt = dxdr[1] * drdt[whichEdge] + dxds[1] * dsdt[whichEdge];
    l += w * sqrt(dxdt * dxdt + dydt * dydt);
  }
  return l;
}


thread_local std::vector<double> GEOCOORD(3 * 6, 0.);
thread_local std::vector<double> POS_VEC(3, 0.);
thread_local std::vector<double> X0(3, 0.);
thread_local std::vector<double> X1(3, 0.);
thread_local std::vector<double> X2(3, 0.);
thread_local std::vector<int> GEO_VERTEX(6);

void feNewRecovery::computeRecoveryAtAllElementDOF(PPR recoveredField, const int index,
                                                   const int iElm,
                                                   std::vector<double> &recoveryAtDOFS,
                                                   bool computeAtP1VerticesOnly)
{
  int nFunctions = _intSpace->getNumFunctions() / _intSpace->getNumComponents();

  // Only fill the first 3 entries (even if triangle is P2)
  // To use when the P2 edge vertices are optimized and moved
  if(computeAtP1VerticesOnly) {
    for(int iDOF = 0; iDOF < 3; ++iDOF) {
      int vertex = _cnc->getVertexConnectivity(iElm, iDOF);
      if(recoveredField == PPR::RECOVERY)
        recoveryAtDOFS[iDOF] = recoveryIndependantTerm[_numTotalRecoveries * vertex + index];
      else
        recoveryAtDOFS[iDOF] = derivativeIndependantTerm[_numTotalDerivatives * vertex + index];
    }
    return;
  }

  int geometryDegree = getGeometricInterpolantDegree(_cnc->getInterpolant());
  int interpolantDegree = _intSpace->getPolynomialDegree();
  if(geometryDegree == interpolantDegree) {
    // Polynomials for the geometry and interpolants match: fast case, just get recovery at (0,0)
    for(int iDOF = 0; iDOF < nFunctions; ++iDOF) {
      int vertex = _cnc->getVertexConnectivity(iElm, iDOF);
      if(recoveredField == PPR::RECOVERY) {
        recoveryAtDOFS[iDOF] = recoveryIndependantTerm[_numTotalRecoveries * vertex + index];
      } else {
        recoveryAtDOFS[iDOF] = derivativeIndependantTerm[_numTotalDerivatives * vertex + index];
      }
    }
  } else {
  // Evaluate the recoveries at (0,0) for the 3 dof at P1 vertices
  // and average the evaluations at edge/element dofs
    _mesh->getCoord(_cnc, iElm, GEOCOORD);

    const std::vector<double> &refCoord = _intSpace->getLcoor();

    for(int iVert = 0; iVert < _cnc->getNumVerticesPerElem(); ++iVert)
      GEO_VERTEX[iVert] = _cnc->getVertexConnectivity(iElm, iVert);

    // Evaluation of each of the 3 vertex polynomials at (0,0)
    for(int iDOF = 0; iDOF < 3; ++iDOF) {
      int vertex = _cnc->getVertexConnectivity(iElm, iDOF);
      if(recoveredField == PPR::RECOVERY)
        recoveryAtDOFS[iDOF] = recoveryIndependantTerm[_numTotalRecoveries * vertex + index];
      else
        recoveryAtDOFS[iDOF] = derivativeIndependantTerm[_numTotalDerivatives * vertex + index];
    }

    double r[3], xLoc[2];
    const std::vector<double> &referenceCoordinates = _intSpace->getLcoor();

    for(int iDOF = 3; iDOF < nFunctions; ++iDOF) {
      // Get reference and physical coordinates of dof
      r[0] = referenceCoordinates[3 * iDOF + 0];
      r[1] = referenceCoordinates[3 * iDOF + 1];
      r[2] = referenceCoordinates[3 * iDOF + 2];
      double refCoordAtDOF[3] = {refCoord[iDOF * 3 + 0], refCoord[iDOF * 3 + 1], 0.};
      _geoSpace->interpolateVectorField(GEOCOORD, refCoordAtDOF, POS_VEC);

      dofLocation dloc = _intSpace->getDOFLocation(iDOF);

      if(dloc == dofLocation::EDGE) {
        // Average the evaluations at the DOF of the two polynomials
        // associated to the edge extremities.

        // Determine which edge and 1D parameter t
        // t is the inverse of the linear parametrization
        double t1D, tol = 1e-6;
        int whichEdge = -1;
        if(r[1] < tol) {
          // Parameterization is r = t, s = 0
          whichEdge = 0;
          t1D = r[0];
        } else if(fabs(r[0] + r[1] - 1.) < tol) {
          // Parameterization is r = 1 - t, s = t
          whichEdge = 1;
          t1D = r[1];
        } else {
          // Parameterization is r = 0, s = 1 - t
          whichEdge = 2;
          t1D = 1. - r[1];
        }

        std::vector<double> dxdr(3, 0.), dxds(3, 0.);

        // Previous/next edge vertex
        int indexEdgeVertex[2];
        if(geometryDegree == 1) {
          indexEdgeVertex[0] = whichEdge;
          indexEdgeVertex[1] = (whichEdge + 1) % 3;
        } else {
          if(t1D < 0.5) {
            indexEdgeVertex[0] = whichEdge;
            indexEdgeVertex[1] = whichEdge + 3;
          } else {
            indexEdgeVertex[0] = whichEdge + 3;
            indexEdgeVertex[1] = (whichEdge + 1) % 3;
          }
        }

        // Evaluate the polynomials of the adjacent vertices
        double pFromAdjacent[2], dist[2];
        for(int iE = 0; iE < 2; ++iE) {
          int ve = GEO_VERTEX[indexEdgeVertex[iE]];
          double xe = GEOCOORD[3 * indexEdgeVertex[iE] + 0];
          double ye = GEOCOORD[3 * indexEdgeVertex[iE] + 1];

          // Local coordinates in the frame centered at vertex ve0/ve1
          xLoc[0] = POS_VEC[0] - xe;
          xLoc[1] = POS_VEC[1] - ye;
          pFromAdjacent[iE] = evaluatePolynomial(recoveredField, index, ve, xLoc);

          // Compute distance to edge extremities
          // if(geometryDegree == 1) {
            // Straight mesh thus euclidian distance for weighting
            dist[iE] =
              sqrt((POS_VEC[0] - xe) * (POS_VEC[0] - xe) + (POS_VEC[1] - ye) * (POS_VEC[1] - ye));
          // }
        }

        // Weighting taking into account the curvature
        // It does not change much and it's slow!
        
        // if(geometryDegree > 1) {
        //   // Compute curved distances from t1D to closest vertices

        //   // double l0, l1;
        //   // computeTransformedEdgeLength(GEOCOORD, dxdr, dxds, t1D, whichEdge, l0, l1);
        //   // dist[0] = l0;
        //   // dist[1] = l1;

        //   if(t1D < 0.5) {
        //     double l0 = computeTransformedEdgeLength2(GEOCOORD, dxdr, dxds, 0., t1D, whichEdge);
        //     double l2 = computeTransformedEdgeLength2(GEOCOORD, dxdr, dxds, 0., 0.5, whichEdge);
        //     double l1 = l2-l0;
        //     dist[0] = l0;
        //     dist[1] = l1;
        //   } else {
        //     double l0 = computeTransformedEdgeLength2(GEOCOORD, dxdr, dxds, 0.5, t1D, whichEdge);
        //     double l2 = computeTransformedEdgeLength2(GEOCOORD, dxdr, dxds, 0.5,  1., whichEdge);
        //     double l1 = l2-l0;
        //     dist[0] = l0;
        //     dist[1] = l1;
        //   }
        // }

        // Averaged evaluations at DOF
        double alpha = dist[1] / (dist[0] + dist[1]);
        recoveryAtDOFS[iDOF] = alpha * pFromAdjacent[0] + (1. - alpha) * pFromAdjacent[1];

      } else if(dloc == dofLocation::ELEMENT) {
        // What weights should we use to average the evaluations
        // at DOFs inside the element?
        // Try just interpolating using the geometry's space for now:
        // This is correct for a P1 mesh, because then it's the barycentric coordinates.
        // For P2 mesh, the Lagrange P2 function can be negative, so it's not a good choice for
        // weights... Try Bézier functions?

        int numFunctionsGeometry = _cnc->getFeSpace()->getNumFunctions();
        std::vector<double> evaluationsAtElementDOF(numFunctionsGeometry);
        for(int i = 0; i < numFunctionsGeometry; ++i) {
          int vertex = _cnc->getVertexConnectivity(iElm, i);
          double xgeo = GEOCOORD[3 * i + 0];
          double ygeo = GEOCOORD[3 * i + 1];
          xLoc[0] = POS_VEC[0] - xgeo;
          xLoc[1] = POS_VEC[1] - ygeo;
          evaluationsAtElementDOF[i] = evaluatePolynomial(recoveredField, index, vertex, xLoc);
        }

        // Weight the evaluations at DOF inside element using Lagrange functions for the geometric
        // space (see remark above, ideally it should be always positive weights).
        recoveryAtDOFS[iDOF] = _geoSpace->interpolateField(evaluationsAtElementDOF, r);
      }
    }
  }
}

double feNewRecovery::averageRecoveryEvaluationsAtQuadNode(PPR recoveredField, const int index,
                                                           const int iElm, const int iQuadNode)
{
  _mesh->getCoord(_cnc, iElm, GEOCOORD);
  _geoSpace->interpolateVectorFieldAtQuadNode(GEOCOORD, iQuadNode, POS_VEC);
  double val = 0., xLoc[2];
  if(_reconstructAtHighOrderNodes) {

    // Option 1: average the 6 recoveries and weight using Lagrange functions
    // Not great because we would need some barycentric-like coordinates to average
    std::vector<double> res(_nNodePerElm, 0.);
    for(int i = 0; i < _nNodePerElm; ++i) {
      int vertex = _cnc->getVertexConnectivity(iElm, i);
      xLoc[0] = POS_VEC[0] - GEOCOORD[3 * i + 0];
      xLoc[1] = POS_VEC[1] - GEOCOORD[3 * i + 1];
      res[i] = evaluatePolynomial(recoveredField, index, vertex, xLoc);
    }

    // Interpolate evaluations? Ideally it should be barycentric weighting?
    val = _intSpace->interpolateFieldAtQuadNode(res, iQuadNode);

    // // Option 2: average 3 recoveries and use barycentric coordinates
    // const std::vector<double> barycentricCoord =
    // _intSpace->getBarycentricCoordinatesAtQuadNode();

    // val = 0.;
    // for(int i = 0; i < 3; ++i) {
    //   int vertex = _cnc->getVertexConnectivity(iElm, i);
    //   xLoc[0] = POS_VEC[0] - GEOCOORD[3 * i + 0];
    //   xLoc[1] = POS_VEC[1] - GEOCOORD[3 * i + 1];
    //   val += barycentricCoord[3 * iQuadNode + i] * evaluatePolynomial(recoveredField, index,
    //   vertex, xLoc);
    // }
    return val;

  } else {
    feInfo("averageRecoveryEvaluationsAtQuadNode not ready");
    exit(-1);
  }

  return val;
}

double feNewRecovery::averageRecoveryEvaluations(PPR recoveredField, const int index,
                                                 const int iElm, double *r)
{
  _mesh->getCoord(_cnc, iElm, GEOCOORD);
  _geoSpace->interpolateVectorField(GEOCOORD, r, POS_VEC);
  double val = 0., xLoc[2];
  if(_reconstructAtHighOrderNodes) {
    // Option 1: average 6 recoveries and use Lagrange functions
    std::vector<double> res(_nNodePerElm, 0.);
    for(int i = 0; i < _nNodePerElm; ++i) {
      int vertex = _cnc->getVertexConnectivity(iElm, i);
      xLoc[0] = POS_VEC[0] - GEOCOORD[3 * i + 0];
      xLoc[1] = POS_VEC[1] - GEOCOORD[3 * i + 1];
      res[i] = evaluatePolynomial(recoveredField, index, vertex, xLoc);
    }
    val = _intSpace->interpolateField(res, r);

    // // Option 2: average 3 recoveries and use barycentric coordinates
    // double barycentricCoord[3] = {1. - r[0] - r[1], r[0], r[1]};

    // val = 0.;
    // for(int i = 0; i < 3; ++i) {
    //   int vertex = _cnc->getVertexConnectivity(iElm, i);
    //   xLoc[0] = POS_VEC[0] - GEOCOORD[3 * i + 0];
    //   xLoc[1] = POS_VEC[1] - GEOCOORD[3 * i + 1];
    //   val += barycentricCoord[i] * evaluatePolynomial(recoveredField, index, vertex, xLoc);
    // }
    return val;

  } else {
    feInfo("averageRecoveryEvaluationsAtQuadNode not ready");
    exit(-1);
  }

  return val;
}

// Interpolates the recovered solution or its derivatives at point x.
// Recovered polynomials are evaluated at the element DOF following Zhang & Naga,
// then the values at DOFs are interpolated using interpolation function from the feSpace
// (not the geometric interpolant)
// Note: this is *not* equivalent to evaluating each polynomial at the quad node, then averaging

thread_local std::vector<double> RECOVERY_AT_DOFS(15);

double feNewRecovery::evaluateRecoveryAtQuadNode(PPR recoveredField, const int index,
                                                 const int iElm, const int iQuadNode)
{
#if defined(AVERAGE_EVALUATIONS_FROM_DOFS)
  // Attempt: evaluating each polynomial at the quad node, then averaging
  return averageRecoveryEvaluationsAtQuadNode(recoveredField, index, iElm, iQuadNode);
#else
  // The Zhang and Naga method: weight evaluation at DOFS then use the Lagrange interpolation
  // functions
  computeRecoveryAtAllElementDOF(recoveredField, index, iElm, RECOVERY_AT_DOFS);
  if(_intSpace->getNumComponents() > 1)
    return _intSpace->interpolateVectorFieldComponentAtQuadNode(RECOVERY_AT_DOFS, iQuadNode,
                                                                _componentToRecover);
  else
    return _intSpace->interpolateFieldAtQuadNode(RECOVERY_AT_DOFS, iQuadNode);
#endif
}

// Evaluate recovered field (recovery or derivative)
// Point x has to be located in the mesh: expensive
double feNewRecovery::evaluateRecovery(PPR recoveredField, const int index, const double *x,
                                       bool averageEvaluations)
{
  int elm;
  double UVW[3];
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("Point (%+-1.6e, %+-1.6e) was not found in the mesh.\n", x[0], x[1]);
    return 0.0;
  } else {
    if(averageEvaluations) {
      return averageRecoveryEvaluations(recoveredField, index, elm, UVW);
    } else {
      // The "classic" Zhang and Naga method: interpolate the data at DOFS using the interpolation
      // functions
      computeRecoveryAtAllElementDOF(recoveredField, index, elm, RECOVERY_AT_DOFS);
      if(_intSpace->getNumComponents() > 1)
        return _intSpace->interpolateVectorFieldComponent(RECOVERY_AT_DOFS, _componentToRecover,
                                                          UVW);
      else
        return _intSpace->interpolateField(RECOVERY_AT_DOFS, UVW);
    }
  }
}

// Evaluate recovery and provide a list of elements to search to locate the point
// to speed up the localization
double feNewRecovery::evaluateRecovery(PPR recoveredField, const int index, const double *x,
                                       const std::vector<int> &elementsToSearch)
{
  int elm;
  double UVW[3];
  bool isFound =
    static_cast<feMesh2DP1 *>(_mesh)->locateVertexInElements(_cnc, x, elementsToSearch, elm, UVW);
  if(!isFound) {
    feWarning("Point (%+-1.6e, %+-1.6e) was not found in the mesh.\n", x[0], x[1]);
    exit(-1);
  } else {
    computeRecoveryAtAllElementDOF(recoveredField, index, elm, RECOVERY_AT_DOFS);
    return _intSpace->interpolateField(RECOVERY_AT_DOFS, UVW);
  }
}

// Evaluate recovered field using only the information from P1 vertices
double feNewRecovery::evaluateRecoveryLinear(PPR recoveredField, const int index, const double *x)
{
  int elm;
  double UVW[3];
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    return 0.0;
  } else {
    // Get recovery only at P1 vertices even is triangle is P2
    computeRecoveryAtAllElementDOF(recoveredField, index, elm, RECOVERY_AT_DOFS, true);
    // Interpolate using barycentric coordinates
    double barycentricCoord[3] = {1. - UVW[0] - UVW[1], UVW[0], UVW[1]};
    double res = 0.;
    for(int i = 0; i < 3; ++i) {
      res += barycentricCoord[i] * RECOVERY_AT_DOFS[i];
    }
    return res;
  }
}

feNewRecovery::feNewRecovery(feSpace *space, int indexComponent, feMesh *mesh, feSolution *sol,
                             std::string meshName, std::string metricMeshName,
                             bool reconstructAtHighOrderNodes,
                             bool useOriginalZhangNagaPatchDefinition, bool append,
                             feMetric *metricField, feMetaNumber *numbering,
                             bool skipRecovery)
  : _mesh(mesh), _sol(sol), _numbering(numbering), _componentToRecover(indexComponent),
    _intSpace(space), _cnc(space->getCncGeo()), _nElm(_cnc->getNumElements()),
    _nNodePerElm(_cnc->getNumVerticesPerElem()), _geoSpace(_cnc->getFeSpace()),
    _degSol(space->getPolynomialDegree()), _degRec(_degSol + 1), _dim(mesh->getDim()),
    _reconstructAtHighOrderNodes(reconstructAtHighOrderNodes)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "SOLUTION AND DERIVATIVES RECOVERY:");
  feInfoCond(FE_VERBOSE > 0, "\t\tRecovering field: \"%s\"", _intSpace->getFieldID().c_str());
  feInfoCond(FE_VERBOSE > 0, "\t\tRecovered fields will be written to file: \"%s\"",
             metricMeshName.c_str());

// Still have to modify the integral formulation for vector FE space
#if !defined(ORIGINAL_ZHANG_NAGA)
  if(space->getNumComponents() > 1) {
    feErrorMsg(
      "Must modify the integral formulation of Zhang and Naga for vector-valued FE spaces first.");
    feErrorMsg("Adjust size of the solution and adr vectors and check indexing.");
    exit(-1);
  }
#endif

  if(reconstructAtHighOrderNodes && getGeometricInterpolantDegree(_cnc->getInterpolant()) == 1) {
    feWarning("The option reconstructAtHighOrderNodes is set to true, but the connectivity is P1.");
    feWarning("Set the mesh to 'curved' to reconstruct on high order nodes.");
    exit(-1);
  }

  // If skipRecovery = true, we only initialize the structure
  // but don't compute the actual recoveries.
  // This is when analytic derivatives are used to adapt.
  // Not very pretty, the data should be taken from elsewhere...
  if(!skipRecovery) {

    _patch = new feNewPatch(_cnc, _mesh, reconstructAtHighOrderNodes,
                            useOriginalZhangNagaPatchDefinition, _degSol, metricField);

    // Assign a 1D quadrature rule to compute transformed edges lengths
    feQuadrature rule1d(8, geometryType::LINE);
    _nQuad1d = rule1d.getNumQuadPoints();
    _wQuad1d = rule1d.getWeights();
    _xQuad1d = rule1d.getXPoints();

    this->setDimensions();
    this->setPolynomialExponents();

    tic();
  #if defined(ORIGINAL_ZHANG_NAGA)
    this->computeVertexLeastSquareMatrices2D();
  #else
    this->computeVertexMassMatrices2D();
  #endif
    feInfoCond(FE_VERBOSE > 0, "\t\tComputed inverse mass/least squares matrices at vertices in %f s",
               toc());

    std::filebuf fbIn, fbOut;
    fbIn.open(meshName, std::ios::in);
    std::istream input(&fbIn);
    fbOut.open(metricMeshName, std::ios::out);
    std::ostream output(&fbOut);

    if(_dim > 1) {
      std::string buffer;
      if(append) {
        // Copy .msh file
        while(getline(input, buffer)) {
          output << buffer << std::endl;
        }
        fbIn.close();
      } else {
        // Copy .msh file except for the possible previous NodeData
        while(getline(input, buffer)) {
          if(buffer == "$NodeData") {
            while(buffer != "$EndNodeData") getline(input, buffer);
            getline(input, buffer);
          }
          output << buffer << std::endl;
        }
        fbIn.close();
      }
    }

    int nRecoveredFields = 0;
    int nRecoveredDerivatives = 0;
    int numRecoveries;
    int maxDerivativeDegree = _degSol + 1;
    _numTotalRecoveries = 0;
    _numTotalDerivatives = 0;

    // Need to loop only to maxDerivativeDegree, and not +1 because we take the derivative of the last
    // recovered So recover up to 2 and take derivative for 3rd order derivative for example.
    for(int iDerivative = 0; iDerivative < maxDerivativeDegree; ++iDerivative) {
      // Reconstruct solution (for i = 0) or derivatives (i > 0)
      numRecoveries = pow(_dim, iDerivative);
      _numTotalRecoveries += numRecoveries;
      _numTotalDerivatives += 2. * numRecoveries;

      tic();
  #if defined(ORIGINAL_ZHANG_NAGA)
      this->computeRHSAndSolve_noIntegral(numRecoveries, nRecoveredFields, nRecoveredDerivatives,
                                          iDerivative);
  #else
      this->computeRHSAndSolve(numRecoveries, nRecoveredFields, nRecoveredDerivatives, iDerivative);
  #endif

      for(int iRecovery = 0; iRecovery < numRecoveries; ++iRecovery) {
        computeDerivative(nRecoveredFields, iRecovery, iDerivative, output);
      }
      feInfoCond(FE_VERBOSE > 0, "\t\tComputed %d derivatives of order %d in %f s", 2 * numRecoveries,
                 iDerivative + 1, toc());

      if(iDerivative > 0) {
        nRecoveredDerivatives += numRecoveries;
      }

      nRecoveredFields += numRecoveries;
    }

    // Condense the independant terms in a vector

    // Allocate a big vector with zeros at edge vertices where
    // we did not reconstruct the fields.
    // This is ugly, but it's just to have things work when we reconstruct
    // only at P1 vertices on P2 meshes, to compare.
    // When comparing, the norms should also be modified to loop only
    // over the first 3 vertices, otherwise the error at edge vertices 
    // is complete garbage.
    // Yes, this is very ugly...
    std::vector<int> allVertices = _cnc->getVerticesConnectivity();
    std::sort(allVertices.begin(), allVertices.end());
    auto last = std::unique(allVertices.begin(), allVertices.end());
    allVertices.erase(last, allVertices.end());
    size_t nAllVertices = allVertices.size();
    recoveryIndependantTerm.resize(nAllVertices * _numTotalRecoveries);
    derivativeIndependantTerm.resize(nAllVertices * _numTotalDerivatives);

    const std::vector<int> vertices = _patch->getVertices();
    size_t nVertices = vertices.size();

    for(size_t i = 0; i < nVertices; ++i) {
      int v = vertices[i];
      for(size_t j = 0; j < _numTotalRecoveries; ++j) {
        recoveryIndependantTerm[_numTotalRecoveries * v + j] = recoveryCoeff.at(v)[j][0];
      }
      for(size_t j = 0; j < _numTotalDerivatives; ++j) {
        derivativeIndependantTerm[_numTotalDerivatives * v + j] = derivativeCoeff.at(v)[j][0];
      }
    }

    if(_dim > 1) fbOut.close();

    // Compute coefficients of the homogeneous error polynomials
    this->computeHomogeneousErrorPolynomials();
  }
}