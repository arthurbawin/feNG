#include "feNG.h"
#include "feMetric.h"
#include "feMetricTools.h"
#include "feNewRecovery.h"
#include "feMatrixInterface.h"

// #include "../contrib/Eigen/Dense"

extern int FE_VERBOSE;

// The names of the reconstructed fields
static std::map<std::pair<int, int>, std::string> suffix = {
  {{0, 0}, ""},     {{1, 0}, "dx"},   {{1, 1}, "dy"},   {{2, 0}, "dxx"},  {{2, 1}, "dxy"},
  {{2, 2}, "dyx"},  {{2, 3}, "dyy"},  {{3, 0}, "dxxx"}, {{3, 1}, "dxxy"}, {{3, 2}, "dxyx"},
  {{3, 3}, "dxyy"}, {{3, 4}, "dyxx"}, {{3, 5}, "dyxy"}, {{3, 6}, "dyyx"}, {{3, 7}, "dyyy"},
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

feNewPatch::feNewPatch(const feCncGeo *cnc, feMesh *mesh, bool reconstructAtHighOrderNodes, feMetric *metricField)
{
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getVerticesConnectivity();

  int nElm = cnc->getNumElements();
  _nVertices = _vertices.size();
  _nNodePerElm = cnc->getNumVerticesPerElem();
  _nEdgePerElm = cnc->getNumEdgesPerElem();

  // If mesh if high-order, keep only the P1 vertices.
  // This is because the Zhang-Naga method reconstructs fields at vertices,
  // then extends to edges and inside by evaluating/averaging.
  int elementOrder = getGeometricInterpolantDegree(cnc->getInterpolant());

  // We have to reconstruct at P2 nodes to compute metrics on
  // the P2 backmesh
  std::unordered_map<int, bool> isOnEdge;
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < 3; ++j) {
      isOnEdge.insert({_vertices[_nNodePerElm * i + j], false});
    }
    for(int j = 3; j < _nNodePerElm; ++j) {
      isOnEdge.insert({_vertices[_nNodePerElm * i + j], true});
    }
  }

  if(!reconstructAtHighOrderNodes) {
    if(elementOrder > 1) {
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
  }

  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  const std::vector<int> &connecNodes = cnc->getVerticesConnectivity();

  // std::unordered_map<int, double> hMax;
  // for(int i = 0; i < nElm; ++i) {
  //   for(int j = 0; j < 6; ++j) {
  //     hMax[connecNodes[_nNodePerElm * i + j]] = 0.;
  //   }
  // }

  // std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
  // for(int i = 0; i < nElm; ++i) {
  //   mesh->getCoord(cnc->getTag(), i, geoCoord);
  //   double x0 = geoCoord[0];
  //   double y0 = geoCoord[1];
  //   double x1 = geoCoord[3];
  //   double y1 = geoCoord[4];
  //   double x2 = geoCoord[6];
  //   double y2 = geoCoord[7];
  //   double l0, l1, l2;

  //   for(int j = 0; j < 3; ++j) {

  //     int vertexP1 = connecNodes[_nNodePerElm * i + j];
  //     int incidentEdge[2] = {j, (j + 2) % 3};

  //     // Compute length of the two incident edges
  //     if(metricField) {
  //       // Use metric field to get a better patch
  //       MetricTensor &M = metricField->getMetricAtSequentialTag(vertexP1);
  //       l0 = sqrt(M(0,0) * (x0 - x1) * (x0 - x1) + 2. * M(0,1) * (x0 - x1) * (y0 - y1) + M(1,1) * (y0 - y1) * (y0 - y1));
  //       l1 = sqrt(M(0,0) * (x1 - x2) * (x1 - x2) + 2. * M(0,1) * (x1 - x2) * (y1 - y2) + M(1,1) * (y1 - y2) * (y1 - y2));
  //       l2 = sqrt(M(0,0) * (x2 - x0) * (x2 - x0) + 2. * M(0,1) * (x2 - x0) * (y2 - y0) + M(1,1) * (y2 - y0) * (y2 - y0));
  //     } else {
  //       // Use Euclidian metric (bad for aniso meshes)
  //       l0 = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
  //       l1 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  //       l2 = sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0));
  //     }
  //     double L[3] = {l0, l1, l2};
    
  //     for(int iEdge = 0; iEdge < 2; ++iEdge) {
  //       // Compute longest connected edge
  //       hMax.at(vertexP1) = fmax(hMax.at(vertexP1), L[ incidentEdge[iEdge] ]);
  //     }
  //   }
  // }

  // // Compute longest connected edge at midnodes
  // for(int i = 0; i < nElm; ++i) {
  //   for(int j = 3; j < 6; ++j) {
  //     int e1 = j - 3;
  //     int e2 = (j - 2) % 3;
  //     hMax.at(connecNodes[_nNodePerElm * i + j]) = fmax(
  //       hMax.at(connecNodes[_nNodePerElm * i + e1]), hMax.at(connecNodes[_nNodePerElm * i + e2]));
  //   }
  // }

  // // Debug check
  // if(true) {
  //   for(auto &p : hMax) {
  //     std::string isIso = metricField ? "aniso" : "iso";
  //     std::string posfile = "cercle" + std::to_string(p.first) + isIso + ".pos";
  //     FILE *file = fopen(posfile.data(), "w");
  //     fprintf(file, "View\"cercle%s%d\"{\n", isIso.data(), p.first);
  //     double x[2] = {0., 0.};
  //     x[0] = mesh->getVertex(p.first)->x();
  //     x[1] = mesh->getVertex(p.first)->y();
  //     if(metricField) {
  //       drawSingleEllipse(file, x, metricField->getMetricAtSequentialTag(p.first), 1, 30);
  //     } else {
  //       double r = p.second;
  //       drawCircle(file, x, r, 30);
  //     }
  //     fprintf(file, "};\n"); fclose(file);
  //   }
  // }

  // if(metricField) {
  //   // Patches of elements within hMax
  //   for(auto &p : hMax) {
  //     double xv = mesh->getVertex(p.first)->x();
  //     double yv = mesh->getVertex(p.first)->y();
  //     for(int i = 0; i < nElm; ++i) {
  //       mesh->getCoord(cnc->getTag(), i, geoCoord);
  //       double x0 = geoCoord[0];
  //       double y0 = geoCoord[1];
  //       double x1 = geoCoord[3];
  //       double y1 = geoCoord[4];
  //       double x2 = geoCoord[6];
  //       double y2 = geoCoord[7];

  //       MetricTensor &M = metricField->getMetricAtSequentialTag(p.first);
  //       double l0 = sqrt(M(0,0) * (x0 - xv) * (x0 - xv) + 2. * M(0,1) * (x0 - xv) * (y0 - yv) + M(1,1) * (y0 - yv) * (y0 - yv));
  //       double l1 = sqrt(M(0,0) * (x1 - xv) * (x1 - xv) + 2. * M(0,1) * (x1 - xv) * (y1 - yv) + M(1,1) * (y1 - yv) * (y1 - yv));
  //       double l2 = sqrt(M(0,0) * (x2 - xv) * (x2 - xv) + 2. * M(0,1) * (x2 - xv) * (y2 - yv) + M(1,1) * (y2 - yv) * (y2 - yv));
  //       double L[3] = {l0, l1, l2};

  //       for(int j = 0; j < 3; ++j) {
  //         if(L[j] <= 1.) {
  //           vertToElems[p.first].insert(i);
  //           break;
  //         }
  //       }
  //     }
  //   }
  // } else {
    // Patches of adjacent elements
    for(int i = 0; i < nElm; ++i) {
      for(int j = 0; j < _nNodePerElm; ++j) {
        vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i);
      }
    }
  // }

  // if(true) {
  //   // Draw patches
  //   std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
  //   for(auto &v : vertToElems) {
  //     std::string posfile = "patch" + std::to_string(v.first) + ".pos";
  //     FILE *file = fopen(posfile.data(), "w");
  //     fprintf(file, "View\"patch\"{\n");
  //     for(auto elem : v.second){
  //       mesh->getCoord(cnc->getTag(), elem, geoCoord);
  //       cnc->writeElementToPOS(file, geoCoord, 1.);
  //     }
  //     fprintf(file, "};\n"); fclose(file);
  //   }
  // }

  if(mesh->getDim() == 1) {
    // The 1D patch associated to vertices is the 2 surrounding elements.
    // Increase the size of the patch if there are too few elements
    for(auto &v : vertToElems) {
      if(v.second.size() < 2) {
        std::vector<int> toAdd;
        for(auto e : v.second) {
          for(int j = 0; j < _nNodePerElm; ++j) {
            for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
              toAdd.push_back(e2);
            }
          }
        }
        for(auto e : toAdd) v.second.insert(e);
      }
    }
  } else if(mesh->getDim() == 2) {
    // Increase the size of the patch if there are too few elements
    // Do not do this for edge vertices otherwise patch becomes too big
    for(auto &v : vertToElems) {
      if(v.second.size() <= 2 && !isOnEdge.at(v.first)) {
        // if(v.second.size() <= 2) {
        std::vector<int> toAdd;
        for(auto e : v.second) {
          // for(int j = 0; j < _nNodePerElm; ++j) {
          for(int j = 0; j < 3; ++j) {
            for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
              toAdd.push_back(e2);
              // v.second.insert(e2);
            }
          }
        }
        for(auto e : toAdd) v.second.insert(e);
      }
    }
  }

  bool debug = false;
  if(debug) {
    // Draw patches
    std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
    for(auto &v : vertToElems) {
      std::string isIso = metricField ? "aniso" : "iso";
      std::string posfile = "patch" + std::to_string(v.first) + isIso + ".pos";
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
    if(_cnc->getGeometry() == geometryType::TRI) {
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
    }
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
    for(size_t i = 0; i < vertices.size(); ++i) {
      int v = vertices[i];
      double xv = _mesh->getVertex(v)->x();
      double yv = _mesh->getVertex(v)->y();

      // Get patch of elements
      const std::set<int> &elemPatch = _patch->getPatch(v);

#if defined(TEST_WITH_EIGEN)
      Eigen::MatrixXd massMatrix = Eigen::MatrixXd::Zero(_dimRecovery, _dimRecovery);
#else
      SquareMatrix massMatrix(_dimRecovery);
#endif

      for(auto elem : elemPatch) {
        _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

        // Loop over quad points and increment mass matrix
        for(int k = 0; k < nQuad; ++k) {
          double jac = J[nQuad * elem + k];

          // Get coordinates in the frame centered at vertex v
          _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
          xLoc[0] = x[0] - xv;
          xLoc[1] = x[1] - yv;

          for(int i = 0; i < _dimRecovery; ++i) {
            monomials[i] = myPow(xLoc[0], _expX_recov[i]) * myPow(xLoc[1], _expY_recov[i]);
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
#if defined(TEST_WITH_EIGEN)
        _inverseMassMatrices[v] = massMatrix.inverse();
#else
        _inverseMassSquareMatrices[v] = massMatrix.inverse();
#endif
      }
    }
  }

  // std::cout << _inverseMassMatrices[3] << std::endl;
  // _inverseMassSquareMatrices[3].print();

  // Mass matrices centered on edges
  // const std::vector<const Edge*> &edges = _mesh->getEdges();
  // #if defined(HAVE_OMP)
  // #pragma omp parallel
  // #endif
  // {
  //   std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
  //   std::vector<double> x(3, 0.0);
  //   std::vector<double> xLoc(3, 0.0);
  //   std::vector<double> monomials(_dimRecovery, 0.);
  //   const Edge* e;

  //   #if defined(HAVE_OMP)
  //   #pragma omp for
  //   #endif
  //   for(size_t i = 0; i < edges.size(); ++i) {
  //     e = edges[i];
  //     double x0 = _mesh->getVertexFromGmshNodeTag(e->getTag(0))->x();
  //     double y0 = _mesh->getVertexFromGmshNodeTag(e->getTag(0))->y();
  //     double x1 = _mesh->getVertexFromGmshNodeTag(e->getTag(1))->x();
  //     double y1 = _mesh->getVertexFromGmshNodeTag(e->getTag(1))->y();
  //     double xc = (x0 + x1) / 2.0;
  //     double yc = (y0 + y1) / 2.0;

  //     // Get patch of elements
  //     const std::set<int> &elemPatch = _patch->getEdgePatch(e->getTag());

  //     Eigen::MatrixXd massMatrix = Eigen::MatrixXd::Zero(_dimRecovery, _dimRecovery);

  //     for(auto elem : elemPatch) {

  //       _mesh->getCoord(_intSpace->getCncGeoTag(), elem, geoCoord);

  //       // Loop over quad points and increment mass matrix
  //       for(int k = 0; k < nQuad; ++k) {
  //         double jac = J[nQuad * elem + k];

  //         // Get coordinates in the frame centered at vertex v
  //         _geoSpace->interpolateVectorFieldAtQuadNode(geoCoord, k, x);
  //         xLoc[0] = x[0] - xc;
  //         xLoc[1] = x[1] - yc;

  //         for(int i = 0; i < _dimRecovery; ++i) {
  //           monomials[i] = myPow( xLoc[0], _expX_recov[i]) * myPow( xLoc[1], _expY_recov[i]);
  //         }

  //         for(int i = 0; i < _dimRecovery; ++i) {
  //           for(int j = 0; j < _dimRecovery; ++j) {
  //             massMatrix(i, j) += jac * w[k] * monomials[i] * monomials[j];
  //           }
  //         }
  //       }
  //     }
  //     #pragma omp critical
  //     _inverseMassMatricesEdges2[e->getTag()] = massMatrix.inverse();
  //   }
  // }
}

void feNewRecovery::computeRHSAndSolve(int numRecoveries, int nRecoveredFields,
                                       int nRecoveredDerivatives, int iDerivative)
{
  int nQuad = _geoSpace->getNumQuadPoints();
  const std::vector<double> &w = _geoSpace->getQuadratureWeights();
  const std::vector<double> &J = _cnc->getJacobians();

  std::vector<int> &vertices = _patch->getVertices();

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
#if defined(TEST_WITH_EIGEN)
    Eigen::VectorXd sol;
#endif
    std::vector<double> u(numRecoveries, 0.);

#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(auto v : vertices) {
      double xv = _mesh->getVertex(v)->x();
      double yv = _mesh->getVertex(v)->y();

#if defined(TEST_WITH_EIGEN)
      std::vector<Eigen::VectorXd> RHS(numRecoveries, Eigen::VectorXd::Zero(_dimRecovery));
#else
      std::vector<Vector> RHS;
      RHS.clear();
      for(int i = 0; i < numRecoveries; ++i) RHS.push_back(Vector(_dimRecovery));
#endif

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

          for(int i = 0; i < _dimRecovery; ++i) {
            monomials[i] = myPow(xLoc[0], _expX_recov[i]) * myPow(xLoc[1], _expY_recov[i]);
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
              // Interpolate the evaluation of the derivative at P1 (or higher order) vertices
              // Bad?
              // u[iComp] = 0.0;
              // for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
              //   int vNode = _cnc->getVertexConnectivity(elem, iVert);
              //   std::vector<double> &du = derivativeCoeff[vNode][nRecoveredDerivatives + iComp];
              //   u[iComp] += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              // }

              // Average the evaluations at edge vertices, then interpolate using
              // Lagrange functions (Zhang & Naga)
              u[iComp] =
                evaluateRecoveryAtQuadNode(PPR::DERIVATIVE, nRecoveredDerivatives + iComp, elem, k);
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
#if defined(TEST_WITH_EIGEN)
        sol = _inverseMassMatrices[v] * RHS[iComp];
#else
        Vector sol = _inverseMassSquareMatrices[v] * RHS[iComp];
#endif
        for(int i = 0; i < _dimRecovery; ++i) {
          recoveryVector[i] = sol(i);
        }
#if defined(HAVE_OMP)
#pragma omp critical
#endif
        recoveryCoeff[v][nRecoveredFields + iComp] = recoveryVector;
      }
    }
  } // pragma omp parallel
}

void feNewRecovery::computeDerivative(int nRecoveredFields, int indexRecovery, int iDerivative,
                                      std::ostream &output)
{
  std::vector<int> &vertices = _patch->getVertices();
  double tol = 0.0;

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
    }
  }

  // Print the solution at first pass
  if(iDerivative == 0) {
    output << "$NodeData\n";
    output << "1\n\"" << _intSpace->getFieldID() << "\"\n1\n3000\n3\n0\n1\n"
           << vertices.size() << "\n";
    for(auto v : vertices) {
      // printf("Writing %d at vertex %d\n", v, _mesh->getVertex(cnt)->getTag());
      // output << _mesh->getVertex(cnt++)->getTag() << " " << recoveryCoeff[v][indexRecovery][0] <<
      // std::endl;
      output << _mesh->getVertex(v)->getTag() << " " << recoveryCoeff[v][indexRecovery][0]
             << std::endl;
      // output << v << " " << recoveryCoeff[v][indexRecovery][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }
  // Also write the derivatives to the common file
  for(int i = 0; i < _dim; ++i) {
    output << "$NodeData\n";
    std::string fieldName = "d" + std::to_string(iDerivative + 1) + _intSpace->getFieldID() +
                            suffix[{iDerivative + 1, _dim * indexRecovery + i}];
    output << "1\n\"" << fieldName << "\"\n1\n3000\n3\n0\n1\n" << vertices.size() << "\n";
    for(auto v : vertices) {
      // output << _mesh->getVertex(cnt++)->getTag() << " " << derivativeCoeff[v][_dim *
      // indexRecovery
      // + i][0] << std::endl;
      output << _mesh->getVertex(v)->getTag() << " "
             << derivativeCoeff[v][_dim * nRecoveredFields + _dim * indexRecovery + i][0]
             << std::endl;
    }
    output << "$EndNodeData\n";
  }
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
      for(int i = 0; i < nDerivatives; ++i) {
        error[target[i]] += derivativeCoeff[v][offset + i][0];
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
             recoveryCoeff[vertex][index][i];
    }
  } else {
    // Evaluate derivative
    for(int i = 0; i < _dimDerivation; ++i) {
      res += myPow(xLoc[0], _expX_deriv[i]) * myPow(xLoc[1], _expY_deriv[i]) *
             derivativeCoeff[vertex][index][i];
    }
  }
  return res;
}

double feNewRecovery::evaluateRecoveryAtVertex(PPR recoveredField, const int index,
                                               const int vertex)
{
  if(recoveredField == PPR::RECOVERY) {
    return recoveryCoeff[vertex][index][0]; // Evaluate recovery
  } else {
    return derivativeCoeff[vertex][index][0]; // Evaluate derivative
  }
}

void feNewRecovery::computeTransformedEdgeLength(std::vector<double> &geoCoord,
                                                 std::vector<double> &dxdr,
                                                 std::vector<double> &dxds, int whichEdge,
                                                 double &firstHalf, double &secondHalf)
{
  double w, t1, t2, l1 = 0., l2 = 0.;
  double R[3] = {0., 0., 0.}, r[3], s[3];
  double drdt[3] = {1., -1., 0.};
  double dsdt[3] = {0., 1., -1.};

  for(size_t i = 0; i < _nQuad1d; ++i) {
    // Map from [-1,1] to parameter space t in [0,1/2] or in [1/2,1]
    t1 = (_xQuad1d[i] + 1.) / 4.;
    t2 = (_xQuad1d[i] + 3.) / 4.;
    w = _wQuad1d[i] / 4.;

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
    l1 += w * sqrt(dxdt * dxdt + dydt * dydt);

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
    l2 += w * sqrt(dxdt * dxdt + dydt * dydt);
  }

  firstHalf = l1;
  secondHalf = l2;
}

void feNewRecovery::computeRecoveryAtAllElementDOF(PPR recoveredField, const int index,
                                                   const int iElm,
                                                   std::vector<double> &recoveryAtDOFS)
{
  std::vector<double> geoCoord(3 * _nNodePerElm, 0.);
  std::vector<double> pos(3, 0.);

  _mesh->getCoord(_cnc, iElm, geoCoord);

  const std::vector<double> &refCoord = _intSpace->getLcoor();
  double refCoordAtDOF[3];

  int nFunctions = _intSpace->getNumFunctions();

  int v0 = _cnc->getVertexConnectivity(iElm, 0);
  int v1 = _cnc->getVertexConnectivity(iElm, 1);
  int v2 = _cnc->getVertexConnectivity(iElm, 2);
  int V[3] = {v0, v1, v2};

  // Evaluation of each of the 3 polynomials at xDOF
  double p[3], shape[3], xLoc[3], val;

  // Coordinates of the P1 vertices
  std::vector<double> x0(3, 0.);
  double xsi0[3] = {0., 0., 0.};
  std::vector<double> x1(3, 0.);
  double xsi1[3] = {1., 0., 0.};
  std::vector<double> x2(3, 0.);
  double xsi2[3] = {0., 1., 0.};
  _geoSpace->interpolateVectorField(geoCoord, xsi0, x0);
  _geoSpace->interpolateVectorField(geoCoord, xsi1, x1);
  _geoSpace->interpolateVectorField(geoCoord, xsi2, x2);

  for(int iDOF = 0; iDOF < nFunctions; ++iDOF) {
    refCoordAtDOF[0] = refCoord[iDOF * 3 + 0];
    refCoordAtDOF[1] = refCoord[iDOF * 3 + 1];
    refCoordAtDOF[2] = refCoord[iDOF * 3 + 2];
    _geoSpace->interpolateVectorField(geoCoord, refCoordAtDOF, pos);

    for(int iv = 0; iv < 3; ++iv) {
      // Local coordinates in the frame centered at vertex iv
      xLoc[0] = pos[0] - geoCoord[3 * iv + 0];
      xLoc[1] = pos[1] - geoCoord[3 * iv + 1];
      p[iv] = evaluatePolynomial(recoveredField, index, V[iv], xLoc);
    }

    if(iDOF < 3) {
      recoveryAtDOFS[iDOF] = p[iDOF];
    } else {
      // Evaluate derivative at edge (or inside, but not implemented) node
      // CANT USE THE GEOMETRIC INTERPOLANT
      // La fonction de Lagrange s'annule au degré de liberté où on veut évaluer
      // Il faut utiliser le ratio de longueur comme proposé par Zhang et Naga
      // _geoSpace->interpolateField(p, 3, refCoordAtDOF, shape, val);

      double l1, l2;
      std::vector<double> dxdr(3, 0.), dxds(3, 0.);

      // Compute distance to edge extremities
      if(iDOF == 3) {
        // Euclidian distance (wrong when elements are curved)
        // l1 = sqrt((pos[0] - x0[0])*(pos[0] - x0[0]) + (pos[1] - x0[1])*(pos[1] - x0[1]));
        // l2 = sqrt((pos[0] - x1[0])*(pos[0] - x1[0]) + (pos[1] - x1[1])*(pos[1] - x1[1]));

        // // Length of transformed reference edge
        computeTransformedEdgeLength(geoCoord, dxdr, dxds, 0, l1, l2);

        double alpha = l2 / (l1 + l2);
        recoveryAtDOFS[iDOF] = alpha * p[0] + (1. - alpha) * p[1];
      }
      if(iDOF == 4) {
        // l1 = sqrt((pos[0] - x1[0])*(pos[0] - x1[0]) + (pos[1] - x1[1])*(pos[1] - x1[1]));
        // l2 = sqrt((pos[0] - x2[0])*(pos[0] - x2[0]) + (pos[1] - x2[1])*(pos[1] - x2[1]));

        // Length of transformed reference edge
        computeTransformedEdgeLength(geoCoord, dxdr, dxds, 1, l1, l2);

        double alpha = l2 / (l1 + l2);
        recoveryAtDOFS[iDOF] = alpha * p[1] + (1. - alpha) * p[2];
      }
      if(iDOF == 5) {
        // l1 = sqrt((pos[0] - x2[0])*(pos[0] - x2[0]) + (pos[1] - x2[1])*(pos[1] - x2[1]));
        // l2 = sqrt((pos[0] - x0[0])*(pos[0] - x0[0]) + (pos[1] - x0[1])*(pos[1] - x0[1]));

        // Length of transformed reference edge
        computeTransformedEdgeLength(geoCoord, dxdr, dxds, 2, l1, l2);

        double alpha = l2 / (l1 + l2);
        recoveryAtDOFS[iDOF] = alpha * p[2] + (1. - alpha) * p[0];
      }
    }
  }
}

// void feNewRecovery::computeRecoveryAtAllElementDOF2(PPR recoveredField, const int index,
//                                                     const int iElm,
//                                                     std::vector<double> &recoveryAtDOFS)
// {
//   int nFunctions = _intSpace->getNumFunctions();
//   double xLoc[3] = {0., 0., 0.};
//   for(int iDOF = 0; iDOF < nFunctions; ++iDOF) {
//     recoveryAtDOFS[iDOF] =
//       evaluatePolynomial(recoveredField, index, _cnc->getVertexConnectivity(iElm, iDOF), xLoc);
//     // recoveryAtDOFS[iDOF] = recoveryCoeff[_cnc->getVertexConnectivity(iElm, iDOF)][index][0];
//   }
// }

void feNewRecovery::computeRecoveryAtAllElementDOF2(PPR recoveredField, const int index,
                                                    const int iElm,
                                                    std::vector<double> &recoveryAtDOFS)
{
  int nFunctions = _intSpace->getNumFunctions();
  for(int iDOF = 0; iDOF < nFunctions; ++iDOF) {
    int vertex = _cnc->getVertexConnectivity(iElm, iDOF);
    if(recoveredField == PPR::RECOVERY)
      recoveryAtDOFS[iDOF] = recoveryCoeff[vertex][index][0];
    else
      recoveryAtDOFS[iDOF] = derivativeCoeff[vertex][index][0];
  }
}

// Interpolates the recovered solution or its derivatives at point x.
// Recovered polynomials are evaluated at the element DOF following Zhang & Naga,
// then the values at DOFs are interpolated using interpolation function from the feSpace
// (not the geometric interpolant)
// Note: this is *not* equivalent to evaluating each polynomial at the quad node, then averaging
double feNewRecovery::evaluateRecoveryAtQuadNode(PPR recoveredField, const int index,
                                                 const int iElm, const int iQuadNode)
{
  std::vector<double> recoveryAtDOFS(_intSpace->getNumFunctions());

  if(_reconstructAtHighOrderNodes)
    computeRecoveryAtAllElementDOF2(recoveredField, index, iElm, recoveryAtDOFS);
  else
    computeRecoveryAtAllElementDOF(recoveredField, index, iElm, recoveryAtDOFS);

  return _intSpace->interpolateFieldAtQuadNode(recoveryAtDOFS, iQuadNode);
}

double feNewRecovery::evaluateRecovery(PPR recoveredField, const int index, const double *x)
{
  int elm;
  double UVW[3];
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    return 0.0;
  } else {
    std::vector<double> recoveryAtDOFS(_intSpace->getNumFunctions());
    if(_reconstructAtHighOrderNodes)
      computeRecoveryAtAllElementDOF2(recoveredField, index, elm, recoveryAtDOFS);
    else
      computeRecoveryAtAllElementDOF(recoveredField, index, elm, recoveryAtDOFS);
    return _intSpace->interpolateField(recoveryAtDOFS, UVW);
  }
}

feNewRecovery::feNewRecovery(feSpace *space, feMesh *mesh, feSolution *sol, std::string meshName,
                             std::string metricMeshName, bool reconstructAtHighOrderNodes,
                             bool append, feMetric *metricField)
  : _mesh(mesh), _sol(sol), _intSpace(space), _cnc(space->getCncGeo()),
    _nElm(_cnc->getNumElements()), _nNodePerElm(_cnc->getNumVerticesPerElem()),
    _geoSpace(_cnc->getFeSpace()), _degSol(space->getPolynomialDegree()), _degRec(_degSol + 1),
    _dim(mesh->getDim()), _reconstructAtHighOrderNodes(reconstructAtHighOrderNodes)
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "SOLUTION AND DERIVATIVES RECOVERY:");
  feInfoCond(FE_VERBOSE > 0, "\t\tRecovered fields will be written to file: \"%s\"",
             metricMeshName.c_str());

  _adr.resize(_intSpace->getNumFunctions());
  _solution.resize(_adr.size());

  _patch = new feNewPatch(_cnc, _mesh, reconstructAtHighOrderNodes, metricField);

  // Assign a 1D quadrature rule to compute transformed edges lengths
  feQuadrature rule1d(8, geometryType::LINE);
  _nQuad1d = rule1d.getNumQuadPoints();
  _wQuad1d = rule1d.getWeights();
  _xQuad1d = rule1d.getXPoints();

  // Set total number of recoveries/derivations
  _nTotalRecoveries = pow(_dim, _degSol); // TODO : A verifier
  _nTotalDerivations = pow(_dim, _degSol + 1);

  this->setDimensions();
  this->setPolynomialExponents();

  tic();
  this->computeVertexMassMatrices2D();
  feInfoCond(FE_VERBOSE > 0, "\t\tComputed inverse mass matrices at vertices in %f s", toc());

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
  for(int iDerivative = 0; iDerivative < _degSol + 1; ++iDerivative) {
    // Reconstruct solution (for i = 0) or derivatives (i > 0)
    numRecoveries = pow(_dim, iDerivative);

    tic();
    this->computeRHSAndSolve(numRecoveries, nRecoveredFields, nRecoveredDerivatives, iDerivative);
    if(iDerivative == 0) {
      feInfoCond(FE_VERBOSE > 0, "\t\tRecovered solution in %f s", toc());
    } else {
      feInfoCond(FE_VERBOSE > 0, "\t\tRecovered %d derivatives of order %d in %f s", numRecoveries,
                 iDerivative, toc());
    }

    tic();
    for(int iRecovery = 0; iRecovery < numRecoveries; ++iRecovery) {
      computeDerivative(nRecoveredFields, iRecovery, iDerivative, output);
    }
    feInfoCond(FE_VERBOSE > 0, "\t\tComputed derivatives in %f s", toc());

    if(iDerivative > 0) {
      nRecoveredDerivatives += numRecoveries;
    }

    nRecoveredFields += numRecoveries;
  }

  if(_dim > 1) fbOut.close();

  // Compute coefficients of the homogeneous error polynomials
  this->computeHomogeneousErrorPolynomials();
}