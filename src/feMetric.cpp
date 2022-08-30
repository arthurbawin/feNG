#include "feMetric.h"
#include "fullMatrix.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#if defined(HAVE_SOPLEX)
#include "soplex.h"
using namespace soplex;
#endif

#include "ellipseToolbox.h"
#include "feMetricTools.h"
// #include "eigen3/Eigen/Eigen"
// #include "eigen3/Eigen/Eigenvalues"
// #include "eigen3/unsupported/Eigen/MatrixFunctions"
#include "../contrib/Eigen/Eigen"
#include "../contrib/Eigen/Eigenvalues"
#include "../contrib/unsupported/Eigen/MatrixFunctions"

#include <assert.h>
#include <iostream>
#include <fstream>

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew);

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions)
  : _recovery(recovery), _options(metricOptions)
{
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is performed on the 
// initial mesh, not on the Gmsh subsitute (gmshModel).
void feMetric::metricScaling()
{
  std::vector<double> &w = _recovery->_geoSpace->getQuadratureWeights();
  std::vector<double> &J = _recovery->_cnc->getJacobians();

  int nQuad = w.size();
  double N = (double)_options.nTargetVertices;
  double p = _options.LpNorm;
  double deg = (double)_recovery->getDegreeSolution();
  double dim = (double)_recovery->getDim();
  double exponent = p * (deg + 1.0) / (2.0 * (p * (deg + 1.0) + dim));

  double I = 0.0;
  for(int iElm = 0; iElm < _recovery->_nElm; ++iElm) {
    for(int k = 0; k < nQuad; ++k) {
      // Interpolate det(Q) at quad nodes
      for(int iNode = 0; iNode < _recovery->_nNodePerElm; ++iNode) {
        int v = _recovery->_cnc->getNodeConnectivity(iElm, iNode);
        I += J[nQuad * iElm + k] * w[k] * _recovery->_geoSpace->getFunctionAtQuadNode(iNode, k) *
             pow(_metrics[v].determinant(), exponent);
      }
    }
  }

  std::vector<int> &vertices = _recovery->getVertices();

  double hMin = _options.hMin, lMax = 1.0 / (hMin * hMin);
  double hMax = _options.hMax, lMin = 1.0 / (hMax * hMax);

  fullMatrix<double> V(3, 3);
  fullVector<double> S(3);
  SVector3 v0, v1, v2;
  for(auto v : vertices) {
    SMetric3 M = _metrics[v];
    feInfo("Initial metric :");
    M.print("");
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), -1.0 / (p * (deg + 1.0) + dim));

    M *= factor;

    M.eig(V, S, false);
    v0 = SVector3(V(0, 0), V(0, 1), V(0, 2));
    v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
    v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));
    M = SMetric3(fmin(lMax, fmax(lMin, S(0))), fmin(lMax, fmax(lMin, S(1))),
                 fmin(lMax, fmax(lMin, S(2))), v0, v1, v2);

    _metrics[v] = M;
    feInfo("Computed scaled metric with factor %1.16e :", factor);
    M.print("");
  }
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is 
// performed on the Gmsh subsitute (gmshModel).
void feMetric::metricScalingFromGmshSubstitute(){
#if defined(HAVE_GMSH)
  double N = (double)_options.nTargetVertices;
  double p = _options.LpNorm;
  double deg = (double)_recovery->getDegreeSolution();
  double dim = (double)_recovery->getDim();
  double exponent = p * (deg + 1.0) / (2.0 * (p * (deg + 1.0) + dim));

  feInfo("%f - %f - %f - %f - %f", N, p, deg, dim, exponent);

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  // Get the jacobians
  std::vector<double> jac, det, pts;
  gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, pts);

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  // Get the nodes
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int dimEntities = -1;
  int tag = -1;
  int includeBoundary = false;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dimEntities, tag, includeBoundary, false);

  // Compute integral of det^exponent
  double I = 0.0;
  double area = 0.0;
  int nQuad = weights.size();

  feInfo("Computing integral of det(Q) - Looping over %d elements tags", elementTags.size());
  feInfo("Computing integral of det(Q) - Looping over %d elements", elementTags[0].size());
  feInfo("det.size = %d", det.size());
  feInfo("jac.size = %d", jac.size());
  feInfo("basisFun.size = %d", basisFunctions.size());
  feInfo("weights.size = %d", nQuad);
  for(size_t iElm = 0; iElm < elementTags[0].size(); iElm++) {
    for(size_t i = 0; i < nQuad; i++) {
      double interpolatedDet = 0.;
      // Interpolate det(Q) at quad nodes
      for(size_t j = 0; j < 3; j++) {
        interpolatedDet += basisFunctions[3*i+j] * pow(_metricsOnGmshModel[elemNodeTags[0][3*iElm+j]].determinant(), exponent);
      }
      I += weights[i] * det[iElm*nQuad + i] * interpolatedDet;
      area += weights[i] * det[iElm*nQuad + i];
    }
  }
  feInfo("Computed integral I = %1.5e", I);
  feInfo("Computed area = %1.5e", area);

  for(size_t i = 0; i < nodeTags.size(); i++){
    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), -1.0 / (p * (deg + 1.0) + dim));

    M *= factor;
    M(2,2) = 1.0;

    _metricsOnGmshModel[nodeTags[i]] = M;
  }
#endif
}

void feMetric::writeSizeFieldSol2D(std::string solFileName)
{
  // Write the size field to a .sol file
  int dim = _recovery->getDim();
  std::vector<int> &vertices = _recovery->getVertices();

  FILE *myfile = fopen(solFileName.c_str(), "w");
  fprintf(myfile, "MeshVersionFormatted 2\n\n");
  fprintf(myfile, "Dimension 2\n\n");
  fprintf(myfile, "SolAtVertices\n");
  fprintf(myfile, "%ld\n", _recovery->getVertices().size());
  fprintf(myfile, "1 3\n\n");

  for(auto v : vertices) {
    SMetric3 M = _metrics[v];
    if(dim == 2)
      fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f\n", M(0, 0), M(0, 1), M(1, 1));
  }

  fprintf(myfile, "End");
  fclose(myfile);
}

void feMetric::writeSizeFieldSol3D(std::string solFileName)
{
  // Write the size field to a .sol file
  int dim = _recovery->getDim();
  std::vector<int> &vertices = _recovery->getVertices();

  FILE *myfile = fopen(solFileName.c_str(), "w");
  fprintf(myfile, "MeshVersionFormatted 2\n\n");
  fprintf(myfile, "Dimension 3\n\n");
  fprintf(myfile, "SolAtVertices\n");
  fprintf(myfile, "%ld\n", _recovery->getVertices().size());
  fprintf(myfile, "1 3\n\n");

  for(auto v : vertices) {
    SMetric3 M = _metrics[v];
    if(dim == 2)
      fprintf(myfile,
              "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n",
              M(0, 0), M(0, 1), M(1, 1), M(0, 2), M(1, 2), M(2, 2));
  }

  fprintf(myfile, "End");
  fclose(myfile);
}

void feMetric::writeSizeFieldGmsh(std::string meshName, std::string metricMeshName)
{
  // Write the size field as NodeData in a copy of the .msh file.
  // For now, the mesh file is simply copied and the NodeData is written at the end.
  feMesh2DP1 *mesh = dynamic_cast<feMesh2DP1 *>(_recovery->_mesh);

  if(mesh->isMeshFileBinary()) {
    printf("In readGmsh : Error - Only reading ASCII files.\n");
    return;
  }
  // if(mesh->getGmshVersion() != 2.2) {
  //   printf("In readGmsh : Error - Anisotropic size field can only be written to a msh 2.2 file
  //   for "
  //          "MMG compatibility.\n");
  //   return;
  // }

  std::filebuf fbIn, fbOut;
  if(fbIn.open(meshName, std::ios::in)) {
    if(fbOut.open(metricMeshName, std::ios::out)) {
      std::istream input(&fbIn);
      std::ostream output(&fbOut);
      std::string buffer;
      // Copy .msh file except for the possible previous NodeData
      while(getline(input, buffer)) {
        if(buffer == "$NodeData") {
          while(buffer != "$EndNodeData") getline(input, buffer);
          getline(input, buffer);
        }
        output << buffer << std::endl;
      }

      std::vector<int> &vertices = _recovery->getVertices();

      // Append the size field
      output << "$NodeData" << std::endl;
      output << 1 << std::endl; // number of string tag
      output << "\"sizeField:metric\"" << std::endl;
      output << 1 << std::endl; // number of real tag
      output << 0 << std::endl;
      output << 3 << std::endl; // number of integer tag
      output << 0 << std::endl;
      output << 9 << std::endl; // metric type ( 1:scalar, 3:vector, 9:tensor)
      output << vertices.size() << std::endl;
      for(auto v : vertices) {
        SMetric3 M = _metrics[v];
        output << v + 1 << " " << M(0, 0) << " " << M(0, 1) << " " << M(0, 2) << " " << M(1, 0)
               << " " << M(1, 1) << " " << M(1, 2) << " " << M(2, 0) << " " << M(2, 1) << " "
               << M(2, 2) << std::endl;
      }
      output << "$EndNodeData" << std::endl;

      fbOut.close();
    }
    fbIn.close();
  }
}

void feMetric::drawEllipsoids(std::string posFile)
{
  FILE *f = fopen(posFile.c_str(), "w");

  fprintf(f, "View \"ellipsesFromSimulation\"{\n");

  double factor = 100.;
  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  std::vector<int> &vertices = _recovery->getVertices();
  for(auto v : vertices) {
    SMetric3 m = _metrics[v];
    Vertex *vv = _recovery->_mesh->getVertex(v);
    getEllipsePoints(factor * m(0, 0), factor * 2.0 * m(0, 1), factor * m(1, 1), vv->x(), vv->y(),
                     xP, yP);
    for(int i = 0; i < nt; ++i) {
      if(i != nt - 1) {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0.,
                xP[i + 1], yP[i + 1], 0., 1, 1);
      } else {
        fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[i], yP[i], 0., xP[0],
                yP[0], 0., 1, 1);
      }
    }
  }

  fprintf(f, "};");
  fclose(f);
}

void feMetric::computeMetricsHechtKuate()
{
  int degSol = _recovery->getDegreeSolution();
  int dimRecovery = _recovery->getDimRecovery();
  std::vector<int> &expX = _recovery->getXExponentsRecovery();
  std::vector<int> &expY = _recovery->getYExponentsRecovery();
  std::map<int, std::vector<double> > &e = _recovery->getErrorCoefficients();
  std::vector<int> &vertices = _recovery->getVertices();

  int nPhi = _options.nPhi;
  std::vector<double> phi(nPhi, 0.); // Error curve discretization
  for(int i = 0; i < nPhi; ++i) {
    phi[i] = i * 2. * M_PI / nPhi;
  }

  double *xC = (double *)malloc(sizeof(double) * nPhi);
  double *yC = (double *)malloc(sizeof(double) * nPhi);
  double *xNew = (double *)malloc(sizeof(double) * nPhi);
  double *yNew = (double *)malloc(sizeof(double) * nPhi);

  double err;
  int cnt;
  for(auto v : vertices) {
    double A, B, C;
    // Construct the discretization of the error curve
    for(int iPhi = 0; iPhi < nPhi; ++iPhi) {
      err = 0.0;
      cnt = 0;
      for(int n = 0; n < dimRecovery; ++n) {
        // Select only exponents of the homogeneous poynomial of order deg+1
        if(expX[n] + expY[n] == degSol + 1) {
          err += e[v][cnt++] * pow(cos(phi[iPhi]), expX[n]) * pow(sin(phi[iPhi]), expY[n]);
        }
      }
      err = fmax(fabs(err), 1e-8);
      xC[iPhi] = pow(1.0 / fabs(err), 1.0 / (degSol + 1)) * cos(phi[iPhi]);
      yC[iPhi] = pow(1.0 / fabs(err), 1.0 / (degSol + 1)) * sin(phi[iPhi]);
    }

    // Compute brute-force metric
    metricHechtKuate(nPhi, xC, yC, A, B, C, 1e-5, xNew, yNew);
    SMetric3 M;
    M.set_m11(A);
    M.set_m21(C);
    M.set_m22(B);
    // Take absolute value of eigenvalues
    fullMatrix<double> V(3, 3);
    fullVector<double> S(3);
    M.eig(V, S, false);
    SVector3 v0(V(0, 0), V(0, 1), V(0, 2));
    SVector3 v1(V(1, 0), V(1, 1), V(1, 2));
    SVector3 v2(V(2, 0), V(2, 1), V(2, 2));
    M = SMetric3(fabs(S(0)), fabs(S(1)), fabs(S(2)), v0, v1, v2);
    _metrics[v] = M;
  }

  // Scale the metric field to obtain the required mesh complexity
  metricScaling();

  FILE *debugFile;
  if(_options.debug) {
    debugFile = fopen("ellipses.pos", "w");
    fprintf(debugFile, "View \" ellipses \"{\n");

    double factor = 100.;
    int nt = 30;
    std::vector<double> xP(nt, 0.);
    std::vector<double> yP(nt, 0.);

    for(auto v : vertices) {
      Vertex *vv = _recovery->_mesh->getVertex(v);
      SMetric3 M = _metrics[v];
      getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), vv->x(), vv->y(),
                       xP, yP);
      for(int j = 0; j < nt; ++j) {
        if(j != nt - 1) {
          fprintf(debugFile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                  xP[j + 1], yP[j + 1], 0., 1, 1);
        } else {
          fprintf(debugFile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                  xP[0], yP[0], 0., 1, 1);
        }
      }
    }

    fprintf(debugFile, "};");
    fclose(debugFile);
  }

  free(xNew);
  free(yNew);
  free(xC);
  free(yC);
}

// The metric field is stored in the view _metricViewTag in the
// gmsh model given in metricOptions.gmshModel
void feMetric::computeMetricsWithDirectionField()
{
// #if defined(HAVE_GMSH)
//   if(!_options.isGmshModelReady) {
//     printf("In feMetric : Error - No Gmsh model available to compute a metric field.\n"
//            "In feMetric : Error - Create a Gmsh geometry and mesh it first.\n");
//     return;
//   }

//   gmsh::model::setCurrent(_options.modelForMetric);

//   // Get the nodes from the Gmsh model
//   std::vector<std::size_t> nodeTags;
//   std::vector<double> coord;
//   std::vector<double> parametricCoord;
//   int dim = -1;
//   int tag = -1;
//   int includeBoundary = true;
//   gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, includeBoundary, false);

//   // Create a view which will contain the metric field
//   _metricViewTag = gmsh::view::add(":metric", 17);
//   feInfo("_metricViewTag = %d", _metricViewTag);
//   int recoveryViewTag = gmsh::view::add("recovery", 21);
//   std::vector<std::vector<double> > metricData;
//   std::vector<std::vector<double> > recoveryData;

//   double xpos[2], g00, g01, g11, C, S;

//   FILE *debugFile1, *debugFile2;
//   if(_options.debug) {
//     debugFile1 = fopen("directionsBeforeSmoothing.pos", "w");
//     debugFile2 = fopen("directionsAfterSmoothing.pos", "w");
//     fprintf(debugFile1, "View \" \"{\n");
//     fprintf(debugFile2, "View \" \"{\n");
//   }

//   // Compute and smooth the direction field
//   std::map<size_t, double> COS;
//   std::map<size_t, double> SIN;

//   std::map<size_t, double> L1;
//   std::map<size_t, double> L2;

//   feInfo("%d noeuds", nodeTags.size());

//   for(size_t i = 0; i < nodeTags.size(); i++) {
//     const double x = coord[3 * i + 0];
//     const double y = coord[3 * i + 1];

//     xpos[0] = x;
//     xpos[1] = y;

//     switch(_options.directionFieldFromDerivativesOfOrder) {
//       case 1:
//         computeDirectionFieldFromGradient(xpos, C, S, 1e-8, _recovery, debugFile1);
//         break;
//       case 2:
//         computeDirectionFieldFromHessian(xpos, C, S, 1e-8, _recovery, debugFile1);
//         break;
//       case 3:
//         // computeDirectionFieldFromThirdOrderDerivatives(x, y, C, S, _rec, debugFile1);
//         break;
//       default:
//         printf(
//           "In feMetric : Error - Cannot compute direction field from derivatives of this order.\n");
//     }

//     COS[nodeTags[i]] = C;
//     SIN[nodeTags[i]] = S;
//   }

//   int smoothMaxIter = 100;
//   double smoothTol = 0.8;
//   smoothDirections(COS, SIN, debugFile2, smoothMaxIter, smoothTol);

//   FILE *debugFile3, *debugFile4, *debugFile5, *debugFile6;
//   if(_options.debug) {
//     fprintf(debugFile1, "};");
//     fclose(debugFile1);
//     fprintf(debugFile2, "};");
//     fclose(debugFile2);
//     debugFile3 = fopen("ellipses.pos", "w");
//     fprintf(debugFile3, "View \" ellipses \"{\n");
//     debugFile4 = fopen("scaledDirections.pos", "w");
//     fprintf(debugFile4, "View \" \"{\n");
//     debugFile5 = fopen("afterScaling.pos", "w");
//     fprintf(debugFile5, "View \" \"{\n");
//     debugFile6 = fopen("afterGradation.pos", "w");
//     fprintf(debugFile6, "View \" \"{\n");
//   }

//   double factor = 100.;
//   int nt = 30;
//   std::vector<double> xP(nt, 0.);
//   std::vector<double> yP(nt, 0.);

//   std::vector<SMetric3> metrics(nodeTags.size());

//   // Compute the principal sizes
//   const int deg = _options.polynomialDegree;
//   const double lMin = _options.hMin;
//   const double lMax = _options.hMax;
//   const double eps = _options.eTargetError;

//   // Min and max eigenvalues based on sizes
//   double lambdaMax = 1./(lMin*lMin);
//   double lambdaMin = 1./(lMax*lMax);

//   // Parameters for the log-simplex metric computation
//   Eigen::Matrix2d Q, V, D = Eigen::Matrix2d::Identity();
//   Eigen::EigenSolver<Eigen::Matrix2d> es;
//   Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
//   bool res;
//   int maxIter = 20;
//   double tol = 5e-1; // Tolerance on the difference Q-Qprev
//   int nThetaPerQuadrant = 50;
//   int numIter;

//   // Set generic LP solver data
//   linearProblem myLP;
//   myLP.problem.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_WARNING);
//   myLP.problem.setIntParam(SoPlex::OBJSENSE,  SoPlex::OBJSENSE_MINIMIZE);
//   DSVector dummycol(0);
//   myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));
//   myLP.problem.addColReal(LPCol(0.0, dummycol, infinity, -infinity));
//   myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));

//   myLP.row.setMax(3);
//   myLP.prim.reDim(3);
//   myLP.lprow = LPRow(3);
//   myLP.lprow.setRhs(infinity);
//   myLP.lprowset = LPRowSet(4*nThetaPerQuadrant,3);

//   for(int i = 0; i < 4*nThetaPerQuadrant; ++i){
//     DSVector row(3);
//     row.add(0, 1.0);
//     row.add(1, 1.0);
//     row.add(2, 1.0);
//     myLP.problem.addRowReal(LPRow(1.0, row, infinity));
//   }

//   for(size_t i = 0; i < nodeTags.size(); i++) {
//     const double x = coord[3 * i + 0];
//     const double y = coord[3 * i + 1];

//     xpos[0] = x;
//     xpos[1] = y;

//     double C = COS[nodeTags[i]];
//     double S = SIN[nodeTags[i]];

//     SMetric3 Mkuate;

//     double l0, l1; // The sizes, not the eigenvalues
//     switch(deg) {
//       case 1: {
//         const double dtt0_ = fabs(dtt(xpos, C, S, _recovery));
//         const double dtt1_ = fabs(dtt(xpos, -S, C, _recovery));
//         l0 = fabs(dtt0_) > 1e-14 ? pow(2.0 * eps / dtt0_, 0.5) : lMax;
//         l1 = fabs(dtt1_) > 1e-14 ? pow(2.0 * eps / dtt1_, 0.5) : lMax;
//         break;
//       }
//       case 2: {
//         const double sizeAlongIsoline = fabs(dttt(xpos, C, S, _recovery, 0));
//         const double sizeAlongGradient = fabs(dttt(xpos, C, S, _recovery, 1));

//         l0 = fabs(sizeAlongIsoline) > 1e-10 ? pow(6.0 * eps / sizeAlongIsoline, 0.3333) : lMax;
//         l1 = fabs(sizeAlongGradient) > 1e-10 ? pow(6.0 * eps / sizeAlongGradient, 0.3333) : lMax;

//         // computeWorstMetric(200, 100, eps, xpos, C, S, l0, l1, _recovery, Mkuate, lMin, lMax);

//         res = computeMetricLogSimplexCurved(xpos, C, S, _recovery, Q, maxIter, nThetaPerQuadrant, tol, numIter, myLP);
//         if(res){
//           feInfo("Metric (%+-1.5e - %+-1.5e - %+-1.5e) found in %2d iterations (log-simplex method)", 
//             Q(0,0), Q(1,0), Q(1,1), numIter);

//           // Bound the eigenvalues of Q
//           es.compute(Q, true);
//           ev = es.eigenvalues();
//           D(0,0) = fmin(lambdaMax, fmax(lambdaMin, ev(0).real() ));
//           D(1,1) = fmin(lambdaMax, fmax(lambdaMin, ev(1).real() ));
//           Q = es.eigenvectors().real() * D * es.eigenvectors().transpose().real();

//         } else{
//           feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e). Setting Q = lmax*I.", x, y);
//           Q = 1. / (lMax * lMax) * Eigen::Matrix2d::Identity();
//         }

//         break;
//       }
//       default: {
//         printf("In feMetric : Error - No metric computation scheme for deg = 0 or deg > 2\n");
//         l0 = lMax;
//         l1 = lMax;
//       }
//     }

//     l0 = std::min(l0, lMax);
//     l0 = std::max(l0, lMin);
//     l1 = std::min(l1, lMax);
//     l1 = std::max(l1, lMin);

//     fprintf(debugFile4, "VP(%g,%g,0){%g,%g,0};", x, y, l0 * C, l0 * S);
//     // fprintf(debugFile4, "VP(%g,%g,0){%g,%g,0};", x, y, -l1 * S, l1 * C);

//     // Eigenvalues
//     double h0 = 1. / (l0 * l0);
//     double h1 = 1. / (l1 * l1);

//     double g00 = C * C * h0 + S * S * h1;
//     double g11 = S * S * h0 + C * C * h1;
//     double g01 = S * C * (h1 - h0);

//     // Check det(M) is positive
//     if(g00 * g11 - g01 * g01 < 1e-10) {
//       printf(
//         "In feMetric : Error - Metric determinant d = %12.12e is negative at node (%4.4e, %4.4e)\n",
//         g00 * g11 - g01 * g01, x, y);
//       g11 = g00 = 1. / lMax / lMax;
//       g01 = 0.0;
//     }

//     // std::vector<double> vMetric(9);
//     // std::vector<double> vRecovery(1);

//     // vMetric[0] = g00;
//     // vMetric[1] = -g01;
//     // vMetric[2] = 0;

//     // vMetric[3] = -g01;
//     // vMetric[4] = g11;
//     // vMetric[5] = 0;

//     // vMetric[6] = 0;
//     // vMetric[7] = 0;
//     // vMetric[8] = 1.0; // export f as well.
//     // vRecovery[0] = f(_recovery, x, y);

//     // metricData.push_back(vMetric);
//     // recoveryData.push_back(vRecovery);

//     SMetric3 M;
//     // M.set_m11(g00);
//     // M.set_m21(-g01);
//     // M.set_m22(g11);
//     M.set_m11(Q(0,0));
//     M.set_m21(Q(0,1));
//     M.set_m22(Q(1,1));

//     double tol = 1e-1;
//     if(_options.debug) {
//       factor = 100.;
//       getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
//       for(int j = 0; j < nt; ++j) {
//         if(j != nt - 1) {
//           fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                   0., xP[j + 1], yP[j + 1], 0., 1, 1);
//         } else {
//           fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                   0., xP[0], yP[0], 0., 1, 1);
//         }
//       }

//       // if(fabs(x - 0.158414) < tol && fabs(y - 0.453362) < tol){
//       //   FILE *myf = fopen("singleEllipseBissection.pos", "w");
//       //   fprintf(myf, "View \" \"{\n");
//       //   for(int j = 0; j < nt; ++j) {
//       //     if(j != nt - 1) {
//       //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//       //               0., xP[j + 1], yP[j + 1], 0., 1, 1);
//       //     } else {
//       //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//       //               0., xP[0], yP[0], 0., 1, 1);
//       //     }
//       //   }
//       //   fprintf(myf, "};");
//       //   fclose(myf);
//       //   feInfo("Printing ellipse at %1.16e - %1.16e with sizes %1.16e - %1.16e", x, y, l0, l1);
//       //   // exit(-1);
//       // }
//     }

//     // TEMPORAIRE : essai
//     metrics[i] = M;
//     // metrics[i] = Mkuate;
//     _metricsOnGmshModel[nodeTags[i]] = M;
//     // _metricsOnGmshModel[i] = Mkuate;
//     // if(Mkuate.determinant() <= 1e-10){
//     //   feInfo("det = %1.16e", Mkuate.determinant());
//     //   Mkuate.print("");
//     //   exit(-1);
//     // }

//     // Pas besoin d'identifier la métrique à un noeud du maillage initial (pas possible même)

//     // double tol = 1e-5;
//     // bool wasFound = false;
//     // std::vector<int> &vertices = _recovery->getVertices();
//     // for(auto v : vertices){
//     //   Vertex *vv = _recovery->_mesh->getVertex(v);
//     //   if(fabs(vv->x() - x) < tol && fabs(vv->y() - y) < tol){
//     //     _metrics[v] = M;
//     //     wasFound = true;
//     //     break;
//     //   }
//     // }
//     // if(!wasFound){
//     //   feInfo("Point %f - %f was not found in step 1", x, y);      
//     //   exit(-1);
//     // }

//     // auto pp = alreadyWrittenTags.insert(nodeTags[i]);

//     // if(pp.second) {
//     //   fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n", nodeTags[i],
//     //   x, y, M(0, 0), M(0, 1), M(1, 1));
//     //   // fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n",
//     //   nodeTags[i],
//     //   // x, y, M(0, 0), M(0, 1), M(1, 1)); Write .sol size field
//     //   // fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n", x, y,
//     //   M(0, 0), M(0, 1), M(1, 1)); fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f\n", M(0,
//     //   0), M(0, 1), M(1, 1));
//     // }
//   }

//   if(_options.nTargetVertices > 0){
//     metricScalingFromGmshSubstitute();
//     // metricScaling();
//   }

//   // Plot the metrics after scaling but before gradation
//   for(size_t i = 0; i < nodeTags.size(); i++) {
//     const double x = coord[3 * i + 0];
//     const double y = coord[3 * i + 1];

//     SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
//     factor = 100.;
//     getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
//     for(int j = 0; j < nt; ++j) {
//       if(j != nt - 1) {
//         fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                 0., xP[j + 1], yP[j + 1], 0., 1, 1);
//       } else {
//         fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                 0., xP[0], yP[0], 0., 1, 1);
//       }
//     }
//   }

//   gradationMetriques(1.6, 3, coord, _metricsOnGmshModel);
//   // feInfo("%d noeuds", nodeTags.size());

//   for(size_t i = 0; i < nodeTags.size(); i++) {
//     const double x = coord[3 * i + 0];
//     const double y = coord[3 * i + 1];

//     SMetric3 M = _metricsOnGmshModel[nodeTags[i]];

//     std::vector<double> vMetric(9);
//     std::vector<double> vRecovery(1);

//     vMetric[0] = M(0,0);
//     vMetric[1] = M(0,1);
//     vMetric[2] = 0;

//     vMetric[3] = M(0,1);
//     vMetric[4] = M(1,1);
//     vMetric[5] = 0;

//     vMetric[6] = 0;
//     vMetric[7] = 0;
//     vMetric[8] = 1.0; // export f as well.
//     vRecovery[0] = f(_recovery, xpos);

//     metricData.push_back(vMetric);
//     recoveryData.push_back(vRecovery);    

//     double tol = 1e-2;
//     if(_options.debug) {
//       factor = 100.;
//       getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
//       for(int j = 0; j < nt; ++j) {
//         if(j != nt - 1) {
//           fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                   0., xP[j + 1], yP[j + 1], 0., 1, 1);
//         } else {
//           fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//                   0., xP[0], yP[0], 0., 1, 1);
//         }
//       }

//       // if(fabs(x - 0.158414) < tol && fabs(y - 0.453362) < tol){
//       //   FILE *myf = fopen("singleEllipseAfterScaling.pos", "w");
//       //   fprintf(myf, "View \" \"{\n");
//       //   for(int j = 0; j < nt; ++j) {
//       //     if(j != nt - 1) {
//       //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//       //               0., xP[j + 1], yP[j + 1], 0., 1, 1);
//       //     } else {
//       //       fprintf(myf, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
//       //               0., xP[0], yP[0], 0., 1, 1);
//       //     }
//       //   }
//       //   fprintf(myf, "};");
//       //   fclose(myf);
//       // }
//     }
//   }

//   if(_options.debug) {
//     fprintf(debugFile3, "};");
//     fclose(debugFile3);
//     fprintf(debugFile4, "};");
//     fclose(debugFile4);
//     fprintf(debugFile5, "};");
//     fclose(debugFile5);
//     fprintf(debugFile6, "};");
//     fclose(debugFile6);
//   }

//   gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags, metricData);
//   gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags, recoveryData);
//   feInfo("Showing the metric and recovery views added to the modelForMetric");

//   std::vector<int> viewTags;
//   gmsh::view::getTags(viewTags);
//   feInfo("There are %d views in the gmsh model : ", viewTags.size());
//   for(auto val : viewTags){
//     feInfo("View %d with index %d", val, gmsh::view::getIndex(val));
//     feInfo("probing view %d", val);
//     std::vector<double> values;
//     double dist;
//     gmsh::view::probe(val, 0.153352, 0.921018, 0.0, values, dist, -1, 9, false, -1.);
//     feInfo("Probed values are : ");
//     for(auto val : values)
//       std::cout<<val<<std::endl;
//   }

//   // gmsh::fltk::run();

//   gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
//   gmsh::view::write(_metricViewTag, _options.metricMeshNameForMMG);
//   gmsh::view::write(recoveryViewTag, _options.recoveryName);

//   // metricOptions.isGmshModelReady = true;

//   // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
//   // gmsh::write("metric.msh");
//   // gmsh::write(metricMeshName);
//   // gmsh::view::write(viewTagF, "sol.msh");

//   // size_t lastindex = metricMeshName.find_last_of(".");
//   // std::string metricMeshNameRoot = metricMeshName.substr(0, lastindex);
//   // std::string toAdapt = metricMeshNameRoot + "_toAdapt.msh";
//   // gmsh::view::write(viewTag, toAdapt);
//   //  gmsh::view::write(viewTagF,"metric.msh",true);

// #else
//   printf("In feMetric : Error - Gmsh is required to create a metric field based on a direction "
//          "field.\n");
//   return;
// #endif
}

void feMetric::computeMetricsLogSimplex()
{
#if defined(HAVE_GMSH) && defined(HAVE_SOPLEX)
  if(!_options.isGmshModelReady) {
    printf("In feMetric : Error - No Gmsh model available to compute a metric field.\n"
           "In feMetric : Error - Create a Gmsh geometry and mesh it first.\n");
    return;
  }

  feInfo("Computing metrics on Gmsh model %s", _options.modelForMetric.c_str());
  gmsh::model::setCurrent(_options.modelForMetric);

  // Get the nodes from the Gmsh model
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int dim = -1;
  int tag = -1;
  int includeBoundary = true;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, includeBoundary, false);

  if(nodeTags.size() == 0){
    feWarning("nodeTags is empty");
    exit(-1);
  }

  std::vector<Vertex> &meshVertices = _recovery->_mesh->getVertices();

  // Create a nodeTags to feVertex->tag map (brute-force for now)

  // On ne peut pas créer le map pour les maillages aniso intermédiaires, puisqu'il n'y a pas de feMesh associé.
  // Il faut recréer le maillage dans la boucle d'intégration pour pouvoir interpoler les métriques.
  // C'est-à-dire : il faut lancer avec -nLoop 0 -nAdapt N, puisque si -nLoop est autre que 0 on boucle sur le maillage aniso
  // avec la même reconstruction.
  double tol = 1e-5;
  _v2n.clear();
  for(Vertex& v : meshVertices){
    _v2n[&v] = -1;
  }
  for(int i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    for(Vertex& v : meshVertices){
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol){
        _v2n[&v] = nodeTags[i];
      }
    }
  }
  // Check
  for(Vertex& v : meshVertices){
    if(_v2n[&v] == -1){
      feInfo("No Gmsh tag was associated to vertex %f - %f", v.x(), v.y());
      exit(-1);
    } 
    // else{
    //   feInfo("Associated tag %d to vertex %f - %f - check : %f - %f ", _v2n[&v], v.x(), v.y(),
    //     coord[3*_v2n[&v]], coord[3*_v2n[&v]+1]);
    // }
  }

  // Create a view which will contain the metric field
  _metricViewTag = gmsh::view::add(":metric", 17);
  feInfo("_metricViewTag = %d", _metricViewTag);
  int recoveryViewTag = gmsh::view::add("recovery", 21);
  std::vector<std::vector<double> > metricData;
  std::vector<std::vector<double> > recoveryData;

  double xpos[2], C, S;

  // Compute and smooth the direction field
  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    xpos[0] = x;
    xpos[1] = y;

    switch(_options.directionFieldFromDerivativesOfOrder) {
      case 1:
        computeDirectionFieldFromGradient(xpos, C, S, 1e-8, _recovery, nullptr);
        break;
      case 2:
        computeDirectionFieldFromHessian(xpos, C, S, 1e-8, _recovery, nullptr);
        break;
      default:
        printf(
          "In feMetric : Error - Cannot compute direction field from derivatives of this order.\n");
    }
    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;
  }

  // int smoothMaxIter = 100;
  // double smoothTol = 0.8;
  // smoothDirections(COS, SIN, nullptr, smoothMaxIter, smoothTol);

  FILE *debugFile3, *debugFile5, *debugFile6;
  if(_options.debug) {
    debugFile3 = fopen("ellipses.pos", "w");
    fprintf(debugFile3, "View \" ellipses \"{\n");
    debugFile5 = fopen("afterScaling.pos", "w");
    fprintf(debugFile5, "View \" \"{\n");
    debugFile6 = fopen("afterGradation.pos", "w");
    fprintf(debugFile6, "View \" \"{\n");
  }

  double factor = 100.;
  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  std::vector<SMetric3> metrics(nodeTags.size());

  // Compute the principal sizes
  const int deg = _options.polynomialDegree;
  const double lMin = _options.hMin;
  const double lMax = _options.hMax;
  const double eps = _options.eTargetError;

  // Min and max eigenvalues based on sizes
  double lambdaMax = 1./(lMin*lMin);
  double lambdaMin = 1./(lMax*lMax);

  // Parameters for the log-simplex metric computation
  linearProblem myLP;
  Eigen::Matrix2d Q, V, D = Eigen::Matrix2d::Identity();
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  bool res;
  int maxIter = 50;
  tol = 1e-1; // Tolerance on the difference Q-Qprev
  int nThetaPerQuadrant = 50;
  int numIter;
  int size = 2 * 4 * nThetaPerQuadrant; // 2 coordonnees * 4 quadrants * nPhi angles/quadrant
  myLP.lvl1.resize(size, 0.);
  myLP.constraints.resize(size, 0.);

  // Set generic LP solver data
  myLP.problem.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
  myLP.problem.setIntParam(SoPlex::OBJSENSE,  SoPlex::OBJSENSE_MINIMIZE);
  DSVector dummycol(0);
  myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));
  myLP.problem.addColReal(LPCol(0.0, dummycol, infinity, -infinity));
  myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));

  myLP.row.setMax(3);
  myLP.prim.reDim(3);
  myLP.lprow = LPRow(3);
  myLP.lprow.setRhs(infinity);
  myLP.lprowset = LPRowSet(4*nThetaPerQuadrant,3);

  for(int i = 0; i < 4*nThetaPerQuadrant; ++i){
    DSVector row(3);
    row.add(0, 1.0);
    row.add(1, 1.0);
    row.add(2, 1.0);
    myLP.problem.addRowReal(LPRow(1.0, row, infinity));
  }

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    xpos[0] = x;
    xpos[1] = y;

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    res = computeMetricLogSimplexCurved(xpos, C, S, _recovery, Q, maxIter, nThetaPerQuadrant, tol, numIter, myLP);
    if(res){
      feInfo("Metric (%+-1.5e - %+-1.5e - %+-1.5e) found in %2d iterations (log-simplex method) (vertex %d/%d) - Target vertices = %d", 
        Q(0,0), Q(1,0), Q(1,1), numIter, i, nodeTags.size(), _options.nTargetVertices);

      // Bound the eigenvalues of Q
      es.compute(Q, true);
      ev = es.eigenvalues();
      D(0,0) = fmin(lambdaMax, fmax(lambdaMin, ev(0).real() ));
      D(1,1) = fmin(lambdaMax, fmax(lambdaMin, ev(1).real() ));
      Q = es.eigenvectors().real() * D * es.eigenvectors().transpose().real();

    } else{
      feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e). Setting Q = lmax*I.", x, y);
      feInfo("C = %f - S = %f", C, S);
      Q = 1. / (lMax * lMax) * Eigen::Matrix2d::Identity();
      exit(-1);
    }

    SMetric3 M;
    M.set_m11(Q(0,0));
    M.set_m21(Q(0,1));
    M.set_m22(Q(1,1));

    metrics[i] = M;
    _metricsOnGmshModel[nodeTags[i]] = M;

    // Plot the metrics raw metrics before scaling
    factor = 100.;
    getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  // Write metrics for Matlab tests
  int NN = 100;
  std::string filename = "metricsForMatlab_a10_n" + std::to_string(NN) + ".txt";
  FILE *fMatlab = fopen(filename.c_str(), "w");
  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double dx = 2./(NN-1);
    double I = round((y + 1.)/dx) + 1;
    I = NN + 1 - I;
    double J = round((x + 1.)/dx) + 1;

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    fprintf(fMatlab, "%d %d %f %f %f %f\n", (int) I, (int) J, M(0,0), M(0,1), M(1,0), M(1,1));
  }
  fclose(fMatlab);

  if(_options.nTargetVertices > 0){
    metricScalingFromGmshSubstitute();
  }

  // Plot the metrics after scaling but before gradation
  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    factor = 100.;
    getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  if(_options.enableGradation){
    gradationMetriques(_options.gradation, 200, coord, _metricsOnGmshModel);
  }

  // Plot the metrics after gradation
  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    factor = 100.;
    getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                0., xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];

    Eigen::Matrix2d m;
    m(0,0) = M(0,0);
    m(0,1) = M(0,1);
    m(1,0) = M(1,0);
    m(1,1) = M(1,1);
    _metricsOnGmshModel_eigen[nodeTags[i]] = m;

    std::vector<double> vMetric(9);
    std::vector<double> vRecovery(1);

    vMetric[0] = M(0,0);
    vMetric[1] = M(0,1);
    vMetric[2] = 0;

    vMetric[3] = M(0,1);
    vMetric[4] = M(1,1);
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.0; // export f as well.
    vRecovery[0] = f(_recovery, xpos);

    metricData.push_back(vMetric);
    recoveryData.push_back(vRecovery);    
  }

  if(_options.debug) {
    fprintf(debugFile3, "};");
    fclose(debugFile3);
    fprintf(debugFile5, "};");
    fclose(debugFile5);
    fprintf(debugFile6, "};");
    fclose(debugFile6);
  }

  gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags, metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags, recoveryData);

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::view::write(_metricViewTag, _options.metricMeshNameForMMG);
  gmsh::view::write(recoveryViewTag, _options.recoveryName);

  // // Write metrics for Matlab tests
  // int NN = 100;
  // std::string filename = "metricsForMatlab_a10_n" + std::to_string(NN) + ".txt";
  // FILE *fMatlab = fopen(filename.c_str(), "w");
  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   const double x = coord[3 * i + 0];
  //   const double y = coord[3 * i + 1];

  //   double dx = 2./(NN-1);
  //   double I = round((y + 1.)/dx) + 1;
  //   I = NN + 1 - I;
  //   double J = round((x + 1.)/dx) + 1;

  //   SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
  //   fprintf(fMatlab, "%d %d %f %f %f %f\n", (int) I, (int) J, M(0,0), M(0,1), M(1,0), M(1,1));
  // }
  // fclose(fMatlab);


  if(false){
  // Test interpolation
  FILE *f = fopen("interpM.pos", "w");
  FILE *f1 = fopen("interpMda1.pos", "w");
  FILE *f2 = fopen("interpMda2.pos", "w");
  FILE *f3 = fopen("interpMonEdge.pos", "w");
  fprintf(f, "View\"interpM\"{\n");
  fprintf(f1, "View\"interpMda1\"{\n");
  fprintf(f2, "View\"interpMda2\"{\n");
  fprintf(f3, "View\"interpMonEdge\"{\n");

  Eigen::Matrix2d myM;
  Eigen::Matrix2d logM11;
  Eigen::Matrix2d M0;
  Eigen::Matrix2d M1;
  Eigen::Matrix2d M2;
  Eigen::Matrix2d sumduda1M;
  Eigen::Matrix2d sumduda2M;
  Eigen::Matrix2d dMda1;
  Eigen::Matrix2d dMda2;


  for(int ielm = 0; ielm < _recovery->_mesh->_elements.size(); ++ielm){
    Triangle *t = _recovery->_mesh->_elements[ielm];
    Vertex *v1 = t->getVertex(0);
    Vertex *v2 = t->getVertex(1);
    Vertex *v3 = t->getVertex(2);

    double a[2] = {0.01, 0.01};
    // Choose an alpha and interpolate derivatives on the edge
    double p0[2] = {v1->x(), v1->y()};
    double p1[2] = {v2->x(), v2->y()};
    double gamma[2] = {p1[0] - p0[0], p1[1] - p0[1]};

    // Interpolate metric at P2 midnode
    double x12[2] = {(p0[0]+p1[0])/2. + a[0], (p0[1]+p1[1])/2. + a[1]};
    // interpolateMetricP1(x12, logM11, M0, M1, M2, sumduda1M, sumduda2M, false);

    // for(int it = 0; it <= 10; ++it){
    //   double tt = (double) it / 10.;
    //   double pos[2] = {p0[0] + gamma[0]*tt + 4.*tt*(1.-tt)*a[0], p0[1] + gamma[1]*tt + 4.*tt*(1.-tt)*a[1]};

    //   interpolateMetricAndDerivativeOnP2Edge(tt, logM11, M0, M1, sumduda1M, sumduda2M, myM, dMda1, dMda2);

    //   factor = 100.;
    //   bool res = getEllipsePoints(factor * myM(0, 0), factor * 2.0 * myM(0, 1), factor * myM(1, 1), pos[0], pos[1], xP, yP);
    //   if(res){
    //     for(int j = 0; j < nt; ++j) {
    //       if(j != nt - 1) {
    //         fprintf(f3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[j + 1], yP[j + 1], 0., 1, 1);
    //       } else {
    //         fprintf(f3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[0], yP[0], 0., 1, 1);
    //       }
    //     }
    //   }
    //   res = getEllipsePoints(factor * dMda1(0, 0), factor * 2.0 * dMda1(0, 1), factor * dMda1(1, 1), pos[0], pos[1], xP, yP);
    //   if(res){
    //     for(int j = 0; j < nt; ++j) {
    //       if(j != nt - 1) {
    //         fprintf(f1, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[j + 1], yP[j + 1], 0., 1, 1);
    //       } else {
    //         fprintf(f1, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[0], yP[0], 0., 1, 1);
    //       }
    //     }
    //   }
    //   res = getEllipsePoints(factor * dMda2(0, 0), factor * 2.0 * dMda2(0, 1), factor * dMda2(1, 1), pos[0], pos[1], xP, yP);
    //   if(res){
    //     for(int j = 0; j < nt; ++j) {
    //       if(j != nt - 1) {
    //         fprintf(f2, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[j + 1], yP[j + 1], 0., 1, 1);
    //       } else {
    //         fprintf(f2, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
    //                 0., xP[0], yP[0], 0., 1, 1);
    //       }
    //     }
    //   }
    // }

    int ni = 3, nj = 3;
    for(int i = 0; i <= ni; ++i){
      for(int j = 0; j <= nj; ++j){
        double rr = (double) i / (double) ni;
        double ss = (double) j / (double) nj;
        if(rr + ss <= 1.){
          double pos[2] = {(1-rr-ss) * v1->x() + rr*v2->x() + ss*v3->x(), 
                           (1-rr-ss) * v1->y() + rr*v2->y() + ss*v3->y()};
          
          interpolateMetricP1(pos, myM, sumduda1M, sumduda2M);

          // myM = myM.exp();
          factor = 100.;
          bool res = getEllipsePoints(factor * myM(0, 0), factor * 2.0 * myM(0, 1), factor * myM(1, 1), pos[0], pos[1], xP, yP);
          if(res){
            for(int j = 0; j < nt; ++j) {
              if(j != nt - 1) {
                fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                        0., xP[j + 1], yP[j + 1], 0., 1, 1);
              } else {
                fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                        0., xP[0], yP[0], 0., 1, 1);
              }
            }
          }

          fprintf(f, "SP(%f,%f,%f){%f};\n", pos[0], pos[1], 0., 1.);
        }
      }
    }
  }

  fprintf(f, "};\n");
  fprintf(f1, "};\n");
  fprintf(f2, "};\n");
  fprintf(f3, "};\n");
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f);
  }
  

#else
  printf("In feMetric : Error - Gmsh is required to create a metric field and SoPlex is required to compute the metric tensors.\n");
  return;
#endif
}

void feMetric::computeMetrics()
{
  switch(_options.computationMethod) {
    case 1:
      computeMetricsHechtKuate();
      break;
    case 2:
      computeMetricsWithDirectionField();
      break;
    case 3 :
    computeMetricsLogSimplex();
      break;
    default:
      printf("In feMetric : Error - Unknown metric computation method.\n");
  }
}

static double UVW[3];
static double GRADUVW[6];
// static Eigen::Matrix2d logM0;
// static Eigen::Matrix2d logM1;
// static Eigen::Matrix2d logM2;

void feMetric::interpolateMetricP1(const double *x, Eigen::Matrix2d &M,
  Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M)
{
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_recovery->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    return;
  } else {

    // Interpolate the log-metric and take exponential
    Triangle *t = _recovery->_mesh->_elements[elm];

    int t0 = _v2n[t->getVertex(0)];
    int t1 = _v2n[t->getVertex(1)];
    int t2 = _v2n[t->getVertex(2)];

    std::map<int, Eigen::Matrix2d>::iterator it = _metricsOnGmshModel_eigen.find(t0);
    if (it == _metricsOnGmshModel_eigen.end()){
      feInfo("Metric not found at %d", t0);
    }
    it = _metricsOnGmshModel_eigen.find(t1);
    if (it == _metricsOnGmshModel_eigen.end()){
      feInfo("Metric not found at %d", t1);
    }
    it = _metricsOnGmshModel_eigen.find(t2);
    if (it == _metricsOnGmshModel_eigen.end()){
      feInfo("Metric not found at %d", t2);
    }

    Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
    Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
    Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

    Eigen::Matrix2d logM0 = M0.log();
    Eigen::Matrix2d logM1 = M1.log();
    Eigen::Matrix2d logM2 = M2.log();

    M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
    M = M.exp();

    /* The barycentric coordinates (u,v,w) are different from UVW above :

        u = 1 - r - s
        v = r
        w = s

        with UVW[0] = r, UVW[1] = s, UVW[2] = t = 0 in 2D
        Unfortunate notation (-:
    */

    // Gradient of the barycentric coordinates w.r.t. P11^l
    double xa = t->getVertex(0)->x();
    double ya = t->getVertex(0)->y();
    double xb = t->getVertex(1)->x();
    double yb = t->getVertex(1)->y();
    double xc = t->getVertex(2)->x();
    double yc = t->getVertex(2)->y();

    // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
    // associated to the barycentric coordinates. Should we use fabs ?
    double A = (xa-xc)*(yb-yc) - (xb-xc)*(ya-yc);

    GRADUVW[0] =  (yb - yc)/A; // grad u = (duda1, duda2)
    GRADUVW[1] = -(xb - xc)/A;
    GRADUVW[2] =  (yc - ya)/A; // grad v
    GRADUVW[3] = -(xc - xa)/A;
    GRADUVW[4] =  (ya - yb)/A; // grad w
    GRADUVW[5] = -(xa - xb)/A;

    sumduda1M = GRADUVW[0] * logM0 + GRADUVW[2] * logM1 + GRADUVW[4] * logM2;
    sumduda2M = GRADUVW[1] * logM0 + GRADUVW[3] * logM1 + GRADUVW[5] * logM2;
  }
}

void feMetric::interpolateMetricP1(const double *x, const double *gammaOrth, Eigen::Matrix2d &M, Eigen::Matrix2d &sumdudaM)
{
  int elm;
  bool isFound = static_cast<feMesh2DP1 *>(_recovery->_mesh)->locateVertex(x, elm, UVW);
  if(!isFound) {
    feWarning("In interpolateMetricP1 : Point (%f, %f) was not found in the mesh.\n", x[0], x[1]);
    M = Eigen::Matrix2d::Identity();
    return;
  } else {

    // Interpolate the log-metric and take exponential
    Triangle *t = _recovery->_mesh->_elements[elm];

    int t0 = _v2n[t->getVertex(0)];
    int t1 = _v2n[t->getVertex(1)];
    int t2 = _v2n[t->getVertex(2)];

    Eigen::Matrix2d M0 = _metricsOnGmshModel_eigen[t0];
    Eigen::Matrix2d M1 = _metricsOnGmshModel_eigen[t1];
    Eigen::Matrix2d M2 = _metricsOnGmshModel_eigen[t2];

    Eigen::Matrix2d logM0 = M0.log();
    Eigen::Matrix2d logM1 = M1.log();
    Eigen::Matrix2d logM2 = M2.log();

    M = (1. - UVW[0] - UVW[1]) * logM0 + UVW[0] * logM1 + UVW[1] * logM2;
    M = M.exp();

    // Gradient of the barycentric coordinates w.r.t. P11^l
    double xa = t->getVertex(0)->x();
    double ya = t->getVertex(0)->y();
    double xb = t->getVertex(1)->x();
    double yb = t->getVertex(1)->y();
    double xc = t->getVertex(2)->x();
    double yc = t->getVertex(2)->y();
    // Area of the triangle, the 1/2 is be cancelled by the area of the subtriangles
    // associated to the barycentric coordinates. Should we use fabs ?
    double A = (xa-xc)*(yb-yc) - (xb-xc)*(ya-yc);
    // Derivative with respect to alpha in 1D, along the bisector
    double gx = gammaOrth[0];
    double gy = gammaOrth[1];
    double duda = (gx*(yb - yc) - gy*(xb - xc)) / A;
    double dvda = (gx*(yc - ya) - gy*(xc - xa)) / A;
    double dwda = (gx*(ya - yb) - gy*(xa - xb)) / A;

    sumdudaM = duda * logM0 + dvda * logM1 + dwda * logM2;
  }
}

void feMetric::interpolateMetricAndDerivativeOnP2Edge(double t,
  const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
  const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
  Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
{
  double phi20 = (1.-t)*(1.-2.*t);
  double phi02 = t*(2.*t - 1.);
  double L = 4. * t * (1.-t); // phi11

  Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

  // lT = T;

  // // Scaling and squaring : find m s.t. ||T||/(2^m) <= 1
  // int m = 0;
  // if(lT.norm() > 1){
  //   m = 1;
  //   while(true){
  //     lT /= pow(2,m);
  //     if(lT.norm() <= 1.)
  //       break;
  //     m++;
  //   }
  //   // feInfo("Found m = %d", m);
  //   // feInfo("norm = %f", lT.norm());
  //   // std::cout<<lT<<std::endl;
  // }

  Eigen::Matrix2d Tn = T;
  M = T.exp();

  Eigen::Matrix2d dTda1 = L * sumduda1M;
  Eigen::Matrix2d dTda2 = L * sumduda2M;
  Eigen::Matrix2d dTnda1 = dTda1;
  Eigen::Matrix2d dTnda2 = dTda2;

  // Compute derivative of exponential by recursion
  // Term for k = 0 is the zero matrix (derivative of I)
  // Term for k = 1 is dTdai
  int kfac = 1;
  dMda1 = dTda1;
  dMda2 = dTda2;

  for(int k = 2; k < 10; ++k){
    dTnda1 = Tn * dTda1 + dTnda1 * T;
    dTnda2 = Tn * dTda2 + dTnda2 * T;

    Tn *= T; // T^n

    kfac *= k;
    dMda1 += 1./((double) kfac) * dTnda1;
    dMda2 += 1./((double) kfac) * dTnda2;
  }

  // // Compute exponential and derivative of lT instead
  // double lambda = pow(2., -m);
  // Tn = lT;
  // M = lT.exp();
  // M = M.pow(1./lambda);

  // dTda1 = L * sumduda1M;
  // dTda2 = L * sumduda2M;
  // dTnda1 = dTda1;
  // dTnda2 = dTda2;

  // // Compute derivative of exponential by recursion
  // // Term for k = 0 is the zero matrix (derivative of I)
  // // Term for k = 1 is dTdai
  // int kfac = 1;
  // dMda1 = dTda1;
  // dMda2 = dTda2;

  // for(int k = 2; k < 10; ++k){
  //   dTnda1 = Tn * dTda1 + dTnda1 * lT;
  //   dTnda2 = Tn * dTda2 + dTnda2 * lT;

  //   Tn *= lT; // T^n

  //   kfac *= k;
  //   dMda1 += 1./(double) kfac * dTnda1;
  //   dMda2 += 1./(double) kfac * dTnda2;
  // }

  // Amat = M.pow(lambda - 1.);
  // Amat = Amat.inverse();
  // Amat *= 1./lambda;

  // dMda1 = Amat * dMda1;
  // dMda2 = Amat * dMda2;
}

void feMetric::interpolateMetricAndDerivativeOnP2Edge(double t,
  const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
  const Eigen::Matrix2d &sumdudaM, Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
{
  double phi20 = (1.-t)*(1.-2.*t);
  double phi02 = t*(2.*t - 1.);
  double L = 4. * t * (1.-t); // phi11

  Eigen::Matrix2d T = phi20 * M20.log() + phi02 * M02.log() + L * M11.log(); // = log M(x(t))

  Eigen::Matrix2d Tn = T;
  M = T.exp();

  Eigen::Matrix2d dTda = L * sumdudaM;
  Eigen::Matrix2d dTnda = dTda;

  // Compute derivative of exponential by recursion
  // Term for k = 0 is the zero matrix (derivative of I)
  // Term for k = 1 is dTdai
  int kfac = 1;
  dMda = dTda;

  for(int k = 2; k < 10; ++k){
    dTnda = Tn * dTda + dTnda * T;
    Tn *= T; // T^n
    kfac *= k;
    dMda += 1./((double) kfac) * dTnda;
  }
}

void interpolateMetricP1Wrapper(void* metric, const double *x, Eigen::Matrix2d &M,
  Eigen::Matrix2d &sumduda1M, Eigen::Matrix2d &sumduda2M)
{
  feMetric *m = static_cast<feMetric*>(metric);
  m->interpolateMetricP1(x, M, sumduda1M, sumduda2M);
}

void interpolateMetricP1Wrapper1D(void* metric, const double *x, const double *gammaOrth, Eigen::Matrix2d &M,
  Eigen::Matrix2d &sumdudaM)
{
  feMetric *m = static_cast<feMetric*>(metric);
  m->interpolateMetricP1(x, gammaOrth, M, sumdudaM);
}

void interpolateMetricAndDerivativeOnP2EdgeWrapper(void* metric, const double t,
  const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
  const Eigen::Matrix2d &sumduda1M, const Eigen::Matrix2d &sumduda2M,
  Eigen::Matrix2d &M, Eigen::Matrix2d &dMda1, Eigen::Matrix2d &dMda2)
{
  feMetric *m = static_cast<feMetric*>(metric);
  m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumduda1M, sumduda2M, M, dMda1, dMda2);
}

void interpolateMetricAndDerivativeOnP2EdgeWrapper1D(void* metric, const double t,
  const Eigen::Matrix2d &M11, const Eigen::Matrix2d &M20, const Eigen::Matrix2d &M02,
  const Eigen::Matrix2d &sumdudaM, Eigen::Matrix2d &M, Eigen::Matrix2d &dMda)
{
  feMetric *m = static_cast<feMetric*>(metric);
  m->interpolateMetricAndDerivativeOnP2Edge(t, M11, M20, M02, sumdudaM, M, dMda);
}