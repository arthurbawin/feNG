#include "feNG.h"
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

feNewPatch::feNewPatch(const feCncGeo *cnc, feMesh *mesh)
{
  // Get unique vertex indices in the node connectivity
  _vertices = cnc->getVerticesConnectivity();
  std::sort(_vertices.begin(), _vertices.end());
  auto last = std::unique(_vertices.begin(), _vertices.end());
  _vertices.erase(last, _vertices.end());

  _nVertices = _vertices.size();
  _nNodePerElm = cnc->getNumVerticesPerElem();
  _nEdgePerElm = cnc->getNumEdgesPerElem();
  const std::vector<int> &connecNodes = cnc->getVerticesConnectivity();

  int nElm = cnc->getNumElements();
  for(int i = 0; i < nElm; ++i) {
    for(int j = 0; j < _nNodePerElm; ++j) {
      vertToElems[connecNodes[_nNodePerElm * i + j]].insert(i); // TODO : A verifier
    }
  }

  switch(mesh->getDim()) {
    case 1:
      // The 1D patch associated to vertices is the 2 surrounding elements.
      // The patch associated to edges is the union of the extremities' patches.
      for(auto e : mesh->_edges) {
        edgeToElems[e.getTag()].insert(vertToElems[e.getTag(0)].begin(),
                                       vertToElems[e.getTag(0)].end());
        edgeToElems[e.getTag()].insert(vertToElems[e.getTag(1)].begin(),
                                       vertToElems[e.getTag(1)].end());
      }
      // Increase the size of the patch if there are too few elements
      for(auto &v : vertToElems) {
        if(v.second.size() < 2) {
          std::cout << "v " << v.first << std::endl;
          std::cout << std::endl;
          for(auto &val : v.second) std::cout << "elem " << val << std::endl;
          std::cout << std::endl;
          std::vector<int> toAdd;
          for(auto val : toAdd) std::cout << "toAdd " << val << std::endl;
          std::cout << std::endl;
          for(auto e : v.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto val : toAdd) std::cout << "toAdd " << val << std::endl;
          std::cout << std::endl;
          for(auto e : toAdd) v.second.insert(e);
        }
      }
      for(auto &v : vertToElems) {
        std::cout << "v " << v.first << std::endl;
        for(auto &val : v.second) std::cout << "elem " << val << std::endl;
        std::cout << std::endl;
      }

      break;
    case 2:

      for(auto e : mesh->_edges) {
        // Insert the patches of both edge vertices to the edge's patch
        int v0 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(0));
        int v1 = mesh->getVertexSequentialTagFromGmshTag(e.getTag(1));
        edgeToElems[e.getTag()].insert(vertToElems[v0].begin(), vertToElems[v0].end());
        edgeToElems[e.getTag()].insert(vertToElems[v1].begin(), vertToElems[v1].end());
      }

      // Increase the size of the patch if there are too few elements
      for(auto &v : vertToElems) {
        if(v.second.size() <= 2) {
          std::vector<int> toAdd;
          for(auto e : v.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto e : toAdd) v.second.insert(e);
        }
      }

      for(auto &p : edgeToElems) {
        if(p.second.size() <= 2) {
          std::vector<int> toAdd;
          for(auto e : p.second) {
            for(int j = 0; j < _nNodePerElm; ++j) {
              for(auto e2 : vertToElems[connecNodes[_nNodePerElm * e + j]]) {
                toAdd.push_back(e2);
                // v.second.insert(e2);
              }
            }
          }
          for(auto e : toAdd) p.second.insert(e);
        }
      }
      break;
    default:
      printf("In feNewPatch : Error - Element patches are only defined for 1D and 2D meshes.\n");
      return;
  }
}

void feNewRecovery::setDimensions()
{
  // The dimension of the polynomial bases for recoveries is the one of degree k+1 :
  if(_dim == 0)
  {
    feErrorMsg(FE_STATUS_ERROR, "No recovery method for connectivity of dimension 0.");
    exit(-1);
  }
  else if(_dim == 1)
  {
    _dimRecovery   = _degRec + 1;
    _dimDerivation = _degSol + 1;
  }
  else if(_dim == 2)
  {
    if(_cnc->getGeometry() == geometryType::TRI)
    {
      _dimRecovery   = (_degRec + 1) * (_degRec + 2) / 2;
      _dimDerivation = (_degSol + 1) * (_degSol + 2) / 2;
    }
    else if(_cnc->getGeometry() == geometryType::QUAD)
    {
      _dimRecovery   = (_degRec + 1) * (_degRec + 1);
      _dimDerivation = (_degSol + 1) * (_degSol + 1);
    }
    else {
      feErrorMsg(FE_STATUS_ERROR, "No polynomial recovery defined for geometry \"%s\".",
             toString(_cnc->getGeometry()).data());
    }
  }
  else if(_dim == 3)
  {
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

  if(_dim == 1)
  {
    // The 1D basis is 1 x x^2 x^3 ...
    for(int i = 0; i < _dimRecovery; ++i)
    {
      _expX_recov[i] = i;
      if(i < _dimDerivation) _expX_deriv[i] = i;
    }
  }
  else if(_dim == 2)
  {
    if(_cnc->getGeometry() == geometryType::TRI)
    {
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
  }
  else {
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
            monomials[i] = myPow( xLoc[0], _expX_recov[i]) * myPow( xLoc[1], _expY_recov[i]);
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

void feNewRecovery::computeRHSAndSolve(int numRecoveries, int nRecoveredFields, int nRecoveredDerivatives, int iDerivative)
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
    std::vector<double> geoCoord(9, 0.);
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
      for(int i = 0; i < numRecoveries; ++i)
        RHS.push_back(Vector(_dimRecovery));
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

          if(iDerivative == 0){
            // Reconstruct the solution : Simply interpolate the solution at quad nodes
            for(int iComp = 0; iComp < numRecoveries; ++iComp) {
              u[iComp] = _intSpace->interpolateFieldAtQuadNode(solution, k);
            }
          }
          else {
            // Reconstruct the derivatives
            // The contributions to the derivative from the vertices must be evaluated and averaged at
            // quad nodes to avoid a trivial solution

            for(int iComp = 0; iComp < numRecoveries; ++iComp) {

              u[iComp] = 0.0;

              // AVERAGE USING GEOMETRIC INTERPOLATION: (not necessarily using edges)
              // ORDER IS LOST FOR P2+

              // Vertices
              for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
                int vNode = _cnc->getVertexConnectivity(elem, iVert);
                std::vector<double> &du = derivativeCoeff[vNode][nRecoveredDerivatives + iComp];
                u[iComp] += _geoSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              }

              // USING REGULAR INTERPOLATION:

              // // Vertices
              // for(int iVert = 0; iVert < _nVertPerElm; ++iVert) {
              //   int vNode = _cnc->getVertexConnectivity(elem, iVert);
              //   std::vector<double> &du = derivativeCoeff[vNode][iComp];
              //   u[iComp] += _intSpace->getFunctionAtQuadNode(iVert, k) * du[0];
              // }
              // // Edges
              // for(int iEdge = 0; iEdge < _nEdgePerElm; ++iEdge) {
              //   int localEdge = fabs(_cnc->getEdgeConnectivity(elem, iEdge));
              //   std::vector<double> &du = derivativeCoeffOnEdges[localEdge][0][iComp];
              //   u[iComp] += _intSpace->getFunctionAtQuadNode(iEdge + 3, k) * du[0];
              // }
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

void feNewRecovery::computeDerivative(int nRecoveredFields, int indexRecovery, int iDerivative, std::ostream &output)
{
  std::vector<int> &vertices = _patch->getVertices();
  double tol = 0.0;

  if(iDerivative == 0) {
    /* Store the independent term of the solution.
    Used in feRecovery::evalDerivative. */
    std::vector<double> solV(vertices.size(), 0.0);
    for(auto v : vertices) solV[v] = recoveryCoeff[v][0][0];
    derivAtVertices.push_back(solV);
  }

  if(_dim == 1) {
    for(auto v : vertices) {
      std::vector<double> &u = recoveryCoeff[v][nRecoveredFields + indexRecovery];
      std::vector<double> dudx(_dimDerivation, 0.);
      for(int i = 0; i < _dimDerivation; ++i) {
        dudx[i] = ((double)i + 1.0) * u[i + 1];
      }
      derivativeCoeff[v][indexRecovery] = dudx;
    }

    // for(auto e : _mesh->_edges) {
    //   std::vector<double> &u = recoveryCoeffOnEdges[e.getTag()][0][indexRecovery];
    //   std::vector<double> dudx(_dimDerivation, 0.);
    //   for(int i = 0; i < _dimDerivation; ++i) {
    //     dudx[i] = ((double)i + 1.0) * u[i + 1];
    //   }
    //   derivativeCoeffOnEdges[e.getTag()][0][indexRecovery] = dudx;
    // }

  } else if(_dim == 2) {
    // To store the derivatives
    std::vector<double> dudxV(vertices.size(), 0.0);
    std::vector<double> dudyV(vertices.size(), 0.0);
    // std::vector<double> dudxE(_mesh->getNumEdges(), 0.0);
    // std::vector<double> dudyE(_mesh->getNumEdges(), 0.0);

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

      dudxV[v] = dudx[0];
      dudyV[v] = dudy[0];
    }

    derivAtVertices.push_back(dudxV);
    derivAtVertices.push_back(dudyV);

    // Exporter les derivees aux sommets (terme indépendant)
    // std::string filename = "derivees_var" + _intSpace->getFieldID() + "_ordre" +
    //                        std::to_string(iDerivative) + "_indRec" + std::to_string(indexRecovery) +
    //                        ".txt";
    // FILE *myfile = fopen(filename.c_str(), "w");
    // // printf("Writing derivatives to %s\n", filename.c_str());
    // for(auto v : vertices) {
    //   fprintf(myfile, "%+-12.12f \t %+-12.12f\n", derivativeCoeff[v][2 * indexRecovery][0],
    //           derivativeCoeff[v][2 * indexRecovery + 1][0]);
    // }
    // fclose(myfile);
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
      // output << _mesh->getVertex(cnt++)->getTag() << " " << derivativeCoeff[v][_dim * indexRecovery
      // + i][0] << std::endl;
      output << _mesh->getVertex(v)->getTag() << " "
             << derivativeCoeff[v][_dim * nRecoveredFields + _dim * indexRecovery + i][0] << std::endl;
    }
    output << "$EndNodeData\n";
  }
}

void numYExponents(int pos, int nLoc, int &ny) {
  if(nLoc == 1){
    return;
  } else{
    if(pos >= nLoc/2) {
      ny++;
      pos -= nLoc/2;
    }
    numYExponents(pos, nLoc/2, ny);
  }
}

void feNewRecovery::computeHomogeneousErrorPolynomials()
{
  if(_dim == 2) {
    std::vector<int> &vertices = _patch->getVertices();
    std::vector<double> error(_degSol + 2, 0.);

    // Number of derivatives of the highest order
    int nDerivatives = pow(2, _degSol+1);
    std::vector<int> target(nDerivatives);

    // High order derivatives contributes to the
    // entry in the error vector according to
    // its y exponent:
    // Example for cubic homogeneous polynomial:
    // error[0] = x^3 * y^0 coefficient
    // error[1] = x^2 * y^1 coefficient
    // error[2] = x^1 * y^2 coefficient
    // error[3] = x^0 * y^3 coefficient
    for(int i = 0; i < nDerivatives; ++i){
      int ny = 0;
      numYExponents(i, nDerivatives, ny);
      target[i] = ny;
    }

    // Number of computed derivatives of lower order
    int offset = pow(2, _degSol+2) - 2 - nDerivatives;

    for(auto v : vertices) {
      std::fill(error.begin(), error.end(), 0.);
      for(int i = 0; i < nDerivatives; ++i){
        error[ target[i] ] += derivativeCoeff[v][offset + i][0];
      }
      errorCoeff[v] = error;
    }
  }
}

double feNewRecovery::evaluatePolynomial(PPR recoveredField, const int index, const int vertex, const double *xLoc)
{
  double res = 0.;
  if(recoveredField == PPR::RECOVERY){
    // Evaluate recovery
    for(int i = 0; i < _dimRecovery; ++i) {
      _basisRecovery[i] = myPow(xLoc[0], _expX_recov[i]) * myPow(xLoc[1], _expY_recov[i]);
      res += _basisRecovery[i] * recoveryCoeff[vertex][index][i];
    }
  } else{
    // Evaluate derivative
    for(int i = 0; i < _dimDerivation; ++i) {
      _basisDerivative[i] = myPow(xLoc[0], _expX_deriv[i]) * myPow(xLoc[1], _expY_deriv[i]);
      res += _basisDerivative[i] * derivativeCoeff[vertex][index][i];
    }
  }
  return res;
}

double feNewRecovery::evaluateRecoveryAtVertex(PPR recoveredField, const int index, const int vertex)
{
  if(recoveredField == PPR::RECOVERY){
    return recoveryCoeff[vertex][index][0]; // Evaluate recovery
  } else{
    return derivativeCoeff[vertex][index][0]; // Evaluate derivative
  }
}

void feNewRecovery::computeRecoveryAtAllElementDOF(PPR recoveredField, const int index, const int iElm)
{
  _mesh->getCoord(_cnc, iElm, _geoCoord);
  
  const std::vector<double> &refCoord = _intSpace->getLcoor();
  double refCoordAtDOF[3];

  int nFunctions = _intSpace->getNumFunctions();

  int v0 = _cnc->getVertexConnectivity(iElm, 0);
  int v1 = _cnc->getVertexConnectivity(iElm, 1);
  int v2 = _cnc->getVertexConnectivity(iElm, 2);
  int V[3] = {v0, v1, v2};

  // Evaluation of each of the 3 polynomials at xDOF
  double p[3], shape[3], xLoc[3], val;

  for(int iDOF = 0; iDOF < nFunctions; ++iDOF)
  { 
    refCoordAtDOF[0] = refCoord[iDOF * 3 + 0];
    refCoordAtDOF[1] = refCoord[iDOF * 3 + 1];
    refCoordAtDOF[2] = refCoord[iDOF * 3 + 2];
    _geoSpace->interpolateVectorField(_geoCoord, refCoordAtDOF, _pos);

    for(int iv = 0; iv < 3; ++iv){
      // Local coordinates in the frame centered at vertex iv
      xLoc[0] = _pos[0] - _geoCoord[3*iv + 0];
      xLoc[1] = _pos[1] - _geoCoord[3*iv + 1];
      p[iv] = evaluatePolynomial(recoveredField, index, V[iv], xLoc);
    }

    // Interpolate using the Geometric interpolant
    _geoSpace->interpolateField(p, 3, refCoordAtDOF, shape, val);
    _FIELD[iDOF] = val;
  }
}

// Interpolates the recovered solution or its derivatives at point x.
double feNewRecovery::evaluateRecoveryAtQuadNode(PPR recoveredField, const int index, const int iElm, const int iQuadNode)
{
  // _mesh->getCoord(_intSpace->getCncGeoTag(), iElm, _geoCoord);
  // // _geoSpace->interpolateVectorFieldAtQuadNode(_geoCoord, iQuadNode, x);
  
  // const std::vector<double> &refCoord = _intSpace->getLcoor();
  // double refCoordAtDOF[3];

  // int nFunctions = _intSpace->getNumFunctions();

  // int v0 = _cnc->getVertexConnectivity(iElm, 0);
  // int v1 = _cnc->getVertexConnectivity(iElm, 1);
  // int v2 = _cnc->getVertexConnectivity(iElm, 2);
  // int V[3] = {v0, v1, v2};

  // // Evaluation of each of the 3 polynomials at xDOF
  // double p[3], shape[3], xLoc[3], val;

  // for(int iDOF = 0; iDOF < nFunctions; ++iDOF)
  // { 
  //   refCoordAtDOF[0] = refCoord[iDOF * 3 + 0];
  //   refCoordAtDOF[1] = refCoord[iDOF * 3 + 1];
  //   refCoordAtDOF[2] = refCoord[iDOF * 3 + 2];
  //   _geoSpace->interpolateVectorField(_geoCoord, refCoordAtDOF, _pos);

  //   for(int iv = 0; iv < 3; ++iv){
  //     // Local coordinates in the frame centered at vertex iv
  //     xLoc[0] = _pos[0] - _geoCoord[3*iv + 0];
  //     xLoc[1] = _pos[1] - _geoCoord[3*iv + 1];
  //     p[iv] = evaluatePolynomial(recoveredField, index, V[iv], xLoc);
  //   }

  //   // Interpolate using the Geometric interpolant
  //   _geoSpace->interpolateField(p, 3, refCoordAtDOF, shape, val);
  //   _FIELD[iDOF] = val;
  // }

  computeRecoveryAtAllElementDOF(recoveredField, index, iElm);
  return _intSpace->interpolateFieldAtQuadNode(_FIELD, iQuadNode);
}

static double UVW[3];

double feNewRecovery::evaluateRecovery(PPR recoveredField, const int index, const double *x)
{
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    return 0.0;
  } else {
    computeRecoveryAtAllElementDOF(recoveredField, index, elm);
    return _intSpace->interpolateField(_FIELD, UVW);
  }
}

feNewRecovery::feNewRecovery(feSpace *space, feMesh *mesh, feSolution *sol,
                       std::string meshName, std::string metricMeshName, bool append)
  : _mesh(mesh)
  , _sol(sol)
  , _intSpace(space)
  , _cnc(space->getCncGeo())
  , _nElm(_cnc->getNumElements())
  , _nNodePerElm(_cnc->getNumVerticesPerElem())
  , _geoSpace(_cnc->getFeSpace())
  , _degSol(space->getPolynomialDegree())
  , _degRec(_degSol + 1)
  , _dim(mesh->getDim())
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "SOLUTION AND DERIVATIVES RECOVERY:");
  feInfoCond(FE_VERBOSE > 0, "\t\tRecovered fields will be written to file: \"%s\"",
      metricMeshName.c_str());

  _adr.resize(_intSpace->getNumFunctions());
  _solution.resize(_adr.size());

  _FIELD.resize(_intSpace->getNumFunctions());
  _geoCoord.resize(3 * _nNodePerElm, 0.);
  _pos.resize(3, 0.);
  
  _patch = new feNewPatch(_cnc, _mesh);

  // Set total number of recoveries/derivations
  _nTotalRecoveries  = pow(_dim, _degSol); // TODO : A verifier
  _nTotalDerivations = pow(_dim, _degSol + 1);

  this->setDimensions();
  this->setPolynomialExponents();

  _basisRecovery.resize(_dimRecovery, 0.);
  _basisDerivative.resize(_dimDerivation, 0.);

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
    if(iDerivative == 0){
      feInfoCond(FE_VERBOSE > 0, "\t\tRecovered solution in %f s", toc());
    } else{
      feInfoCond(FE_VERBOSE > 0, "\t\tRecovered %d derivatives of order %d in %f s", numRecoveries, iDerivative, toc());
    }

    tic();
    for(int iRecovery = 0; iRecovery < numRecoveries; ++iRecovery) {
      computeDerivative(nRecoveredFields, iRecovery, iDerivative, output);
    }
    feInfoCond(FE_VERBOSE > 0, "\t\tComputed derivatives in %f s", toc());

    if(iDerivative > 0){
      nRecoveredDerivatives += numRecoveries;
    }

    nRecoveredFields += numRecoveries;
  }

  if(_dim > 1) fbOut.close();

  // Compute coefficients of the homogeneous error polynomials
  this->computeHomogeneousErrorPolynomials();
}