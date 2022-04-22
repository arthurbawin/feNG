#include "feMetric.h"
#include "fullMatrix.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

#include "ellipseToolbox.h"
#include "feMetricTools.h"

#include <assert.h>
#include <iostream>
#include <fstream>

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew);

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions)
  : _recovery(recovery), _options(metricOptions)
{
}

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

    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), -1.0 / (p * (deg + 1.0) + dim));

    M *= factor;

    M.eig(V, S, false);
    v0 = SVector3(V(0, 0), V(0, 1), V(0, 2));
    v1 = SVector3(V(1, 0), V(1, 1), V(1, 2));
    v2 = SVector3(V(2, 0), V(2, 1), V(2, 2));
    M = SMetric3(fmin(lMax, fmax(lMin, S(0))), fmin(lMax, fmax(lMin, S(1))),
                 fmin(lMax, fmax(lMin, S(2))), v0, v1, v2);

    _metrics[v] = M;
  }
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

  double factor = 1.;
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

    double factor = 1.;
    int nt = 30;
    std::vector<double> xP(nt, 0.);
    std::vector<double> yP(nt, 0.);

    for(auto v : vertices){
      Vertex *vv = _recovery->_mesh->getVertex(v);
      SMetric3 M = _metrics[v];
      getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), vv->x(), vv->y(), xP, yP);
      for(int j = 0; j < nt; ++j) {
        if(j != nt - 1) {
          fprintf(debugFile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                  0., xP[j + 1], yP[j + 1], 0., 1, 1);
        } else {
          fprintf(debugFile, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j],
                  0., xP[0], yP[0], 0., 1, 1);
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
#if defined(HAVE_GMSH)
  if(!_options.isGmshModelReady) {
    printf("In feMetric : Error - No Gmsh model available to compute a metric field.\n"
           "In feMetric : Error - Create a Gmsh geometry and mesh it first.\n");
    return;
  }

  gmsh::initialize();
  gmsh::model::setCurrent(_options.gmshModel);

  // Get the nodes from the Gmsh model
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int dim = 2;
  int tag = -1;
  int includeBoundary = true;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, includeBoundary, false);
  // Create a view which will contain the metric field
  _metricViewTag = gmsh::view::add(":metric");
  int recoveryViewTag = gmsh::view::add("recovery");
  std::vector<std::vector<double> > metricData;
  std::vector<std::vector<double> > recoveryData;

  double g00, g01, g11, C, S;

  FILE *debugFile1, *debugFile2;
  if(_options.debug) {
    debugFile1 = fopen("directionsBeforeSmoothing.pos", "w");
    debugFile2 = fopen("directionsAfterSmoothing.pos", "w");
    fprintf(debugFile1, "View \" \"{\n");
    fprintf(debugFile2, "View \" \"{\n");
  }

  // Compute and smooth the direction field
  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    switch(_options.directionFieldFromDerivativesOfOrder) {
      case 1:
        computeDirectionFieldFromGradient(x, y, C, S, 1e-8, _recovery, debugFile1);
        break;
      case 2:
        computeDirectionFieldFromHessian(x, y, C, S, 1e-8, _recovery, debugFile1);
        break;
      case 3:
        // computeDirectionFieldFromThirdOrderDerivatives(x, y, C, S, _rec, debugFile1);
        break;
      default:
        printf(
          "In feMetric : Error - Cannot compute direction field from derivatives of this order.\n");
    }

    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;
  }

  int smoothMaxIter = 100;
  double smoothTol = 0.1;
  smoothDirections(COS, SIN, debugFile2, smoothMaxIter, smoothTol);

  FILE *debugFile3, *debugFile4;
  if(_options.debug) {
    fprintf(debugFile1, "};");
    fclose(debugFile1);
    fprintf(debugFile2, "};");
    fclose(debugFile2);
    debugFile3 = fopen("ellipses.pos", "w");
    fprintf(debugFile3, "View \" ellipses \"{\n");
    debugFile4 = fopen("scaledDirections.pos", "w");
    fprintf(debugFile4, "View \" \"{\n");
  }

  double factor = 10.;
  int nt = 30;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  std::vector<SMetric3> metrics(nodeTags.size());

  // // Write the .sol size field file
  // FILE *myfile = fopen("toIntersect2D.sol", "w");
  // fprintf(myfile, "MeshVersionFormatted 2\n\n");
  // fprintf(myfile, "Dimension 2\n\n");
  // fprintf(myfile, "SolAtVertices\n");
  // std::set<int> s(nodeTags.begin(), nodeTags.end());
  // fprintf(myfile, "%ld\n", s.size());
  // fprintf(myfile, "1 3\n\n");
  // FILE *scaled = fopen("scaledDirections.pos", "w");
  // fprintf(scaled, "View \" \"{\n");

  // std::set<std::size_t> alreadyWrittenTags;

  // Compute the principal sizes
  const int deg = _options.polynomialDegree;
  const double lMin = _options.hMin;
  const double lMax = _options.hMax;
  const double eps = _options.eTargetError;

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    double l0, l1; // The sizes, not the eigenvalues
    switch(deg) {
      case 1: {
        const double dtt0_ = fabs(dtt(x, y, C, S, _recovery));
        const double dtt1_ = fabs(dtt(x, y, -S, C, _recovery));
        l0 = fabs(dtt0_) > 1e-14 ? pow(2.0 * eps / dtt0_, 0.5) : lMax;
        l1 = fabs(dtt1_) > 1e-14 ? pow(2.0 * eps / dtt1_, 0.5) : lMax;
        break;
      }
      case 2: {
        const double sizeAlongIsoline = fabs(dttt(x, y, C, S, _recovery, 0));
        // const double sizeAlongGradient = fabs(dttt(x, y, -S, C, _recovery, 1));
        const double sizeAlongGradient = fabs(dttt(x, y, C, S, _recovery, 1));
        l0 = fabs(sizeAlongIsoline) > 1e-14 ? pow(6.0 * eps / sizeAlongIsoline, 0.3333) : lMax;
        l1 = fabs(sizeAlongGradient) > 1e-14 ? pow(6.0 * eps / sizeAlongGradient, 0.3333) : lMax;
        break;
      }
      default: {
        printf("In feMetric : Error - No metric computation scheme for deg = 0 or deg > 2\n");
        l0 = lMax;
        l1 = lMax;
      }
    }

    l0 = std::min(l0, lMax);
    l0 = std::max(l0, lMin);
    l1 = std::min(l1, lMax);
    l1 = std::max(l1, lMin);

    fprintf(debugFile4, "VP(%g,%g,0){%g,%g,0};", x, y, l0 * C, l0 * S);
    fprintf(debugFile4, "VP(%g,%g,0){%g,%g,0};", x, y, -l1 * S, l1 * C);

    // Eigenvalues
    double h0 = 1. / (l0 * l0);
    double h1 = 1. / (l1 * l1);

    double g00 = C * C * h0 + S * S * h1;
    double g11 = S * S * h0 + C * C * h1;
    double g01 = S * C * (h1 - h0);

    // Check det(M) is positive
    if(g00 * g11 - g01 * g01 < 1e-10) {
      printf(
        "In feMetric : Error - Metric determinant d = %12.12e is negative at node (%4.4e, %4.4e)\n",
        g00 * g11 - g01 * g01, x, y);
      g11 = g00 = 1. / lMax / lMax;
      g01 = 0.0;
    }

    std::vector<double> vMetric(9);
    std::vector<double> vRecovery(1);

    vMetric[0] = g00;
    vMetric[1] = -g01;
    vMetric[2] = 0;

    vMetric[3] = -g01;
    vMetric[4] = g11;
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.0; // export f as well.
    vRecovery[0] = f(_recovery, x, y);

    metricData.push_back(vMetric);
    recoveryData.push_back(vRecovery);

    SMetric3 M;
    M.set_m11(g00);
    M.set_m21(-g01);
    M.set_m22(g11);
    metrics[i] = M;

    double tol = 1e-5;
    std::vector<int> &vertices = _recovery->getVertices();
    for(auto v : vertices){
      Vertex *vv = _recovery->_mesh->getVertex(v);
      if(fabs(vv->x() - x) < tol && fabs(vv->y() - y) < tol){
        _metrics[v] = M;
      }
    }   

    if(_options.debug) {
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

    // auto pp = alreadyWrittenTags.insert(nodeTags[i]);

    // if(pp.second) {
    //   fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n", nodeTags[i],
    //   x, y, M(0, 0), M(0, 1), M(1, 1));
    //   // fprintf(ffff, "%d \t %+-4.4e \t %+-4.4e : %+-4.4e \t %+-4.4e \t %+-4.4e \n",
    //   nodeTags[i],
    //   // x, y, M(0, 0), M(0, 1), M(1, 1)); Write .sol size field
    //   // fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f \t %+-10.12f\n", x, y,
    //   M(0, 0), M(0, 1), M(1, 1)); fprintf(myfile, "%+-10.12f \t %+-10.12f \t %+-10.12f\n", M(0,
    //   0), M(0, 1), M(1, 1));
    // }
  }

  if(_options.debug) {
    fprintf(debugFile3, "};");
    fclose(debugFile3);
    fprintf(debugFile4, "};");
    fclose(debugFile4);
  }

  gmsh::view::addModelData(_metricViewTag, 0, _options.gmshModel, "NodeData", nodeTags, metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.gmshModel, "NodeData", nodeTags,
                           recoveryData);
  // gmsh::fltk::run();

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::view::write(_metricViewTag, "metric.msh");
  gmsh::view::write(recoveryViewTag, "recovery.msh");

  // metricOptions.isGmshModelReady = true;

  // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  // gmsh::write("metric.msh");
  // gmsh::write(metricMeshName);
  // gmsh::view::write(viewTagF, "sol.msh");

  // size_t lastindex = metricMeshName.find_last_of(".");
  // std::string metricMeshNameRoot = metricMeshName.substr(0, lastindex);
  // std::string toAdapt = metricMeshNameRoot + "_toAdapt.msh";
  // gmsh::view::write(viewTag, toAdapt);
  //  gmsh::view::write(viewTagF,"metric.msh",true);

#else
  printf("In feMetric : Error - Gmsh is required to create a metric field based on a direction "
         "field.\n");
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
    default:
      printf("In feMetric : Error - Unknown metric computation method.\n");
  }
}