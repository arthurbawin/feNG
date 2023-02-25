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
#include "../contrib/Eigen/Eigen"
#include "../contrib/Eigen/Eigenvalues"

void metricHechtKuate(int nbpoints, double *x, double *y, double &A, double &B, double &C,
                      double epsilon, double *xNew, double *yNew);

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions)
  : _recovery(recovery), _options(metricOptions)
{
}

feMetric::feMetric(feNewRecovery *recovery, feMetricOptions metricOptions)
  : _newRecovery(recovery), _options(metricOptions)
{
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is performed on the
// initial mesh, not on the Gmsh subsitute (gmshModel).
void feMetric::metricScaling()
{
  std::vector<double> &w = _recovery->_geoSpace->getQuadratureWeights();
  const std::vector<double> &J = _recovery->_cnc->getJacobians();

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
        int v = _recovery->_cnc->getVertexConnectivity(iElm, iNode);
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
// template <class MetricType>
void feMetric::metricScalingFromGmshSubstitute()
{
#if defined(HAVE_GMSH)
  double N = (double) _options.nTargetVertices;
  double p = _options.LpNorm;
  double deg = (double) _recovery->getDegreeSolution();
  double dim = (double) _recovery->getDim();

  // double exponent = p * (deg + 1.0) / (2.0 * (p * (deg + 1.0) + dim));
  double exponent = p * (deg + 1.0) / (p * (deg + 1.0) + dim);

  feInfo("%f - %f - %f - %f - %f", N, p, deg, dim, exponent);

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents, basisFunctions,
                                       numOrientations);

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
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dimEntities, tag, includeBoundary,
                              false);

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
        interpolatedDet +=
          basisFunctions[3 * i + j] *
          pow(_metricsOnGmshModel[elemNodeTags[0][3 * iElm + j]].determinant(), exponent);
      }
      I += weights[i] * det[iElm * nQuad + i] * interpolatedDet;
      area += weights[i] * det[iElm * nQuad + i];
    }
  }
  feInfo("Computed integral I = %1.5e", I);
  feInfo("Computed area = %1.5e", area);

  for(size_t i = 0; i < nodeTags.size(); i++) {
    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), -1.0 / (p * (deg + 1.0) + dim));

    M *= factor;
    M(2, 2) = 1.0;

    _metricsOnGmshModel[nodeTags[i]] = M;
  }
#endif
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is
// performed on the Gmsh subsitute (gmshModel).
template <class MetricType>
void feMetric::metricScalingFromGmshSubstitute(std::map<int, MetricType> &metrics,
    const std::vector<size_t> &nodeTags, 
    double exponentInIntegral, double exponentForDeterminant)
{
#if defined(HAVE_GMSH)

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(triP1, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents, basisFunctions,
                                       numOrientations);

  // Get the jacobians
  std::vector<double> jac, det, pts;
  gmsh::model::mesh::getJacobians(triP1, localCoord, jac, det, pts);

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  // Compute integral of det^exponent
  double I = 0.0;
  double area = 0.0;
  int nQuad = weights.size();

  double N = (double) _options.nTargetVertices;
  double dim = (double) _newRecovery->getDim();

  feInfo("Computing integral of det(Q) - Looping over %d elements tags", elementTags.size());
  feInfo("Computing integral of det(Q) - Looping over %d elements", elementTags[0].size());
  feInfo("There are %d quad nodes per elem", nQuad);
  feInfo("det.size = %d", det.size());
  feInfo("jac.size = %d", jac.size());
  feInfo("basisFun.size = %d", basisFunctions.size());
  feInfo("weights.size = %d", nQuad);
  feInfo("Target = %f", N);

  for(size_t iElm = 0; iElm < elementTags[0].size(); iElm++) {

    for(size_t i = 0; i < nQuad; i++) {

      double interpolatedDet = 0.;
      // Interpolate det(H) or det(Q) at quad nodes
      for(size_t j = 0; j < 3; j++) {

        // Get vertex associated to Gmsh local nodetag
        // const Vertex *v = _n2v[ elemNodeTags[0][3 * iElm + j] ];
        // interpolatedDet += basisFunctions[3 * i + j]
        //  * pow( metrics[v].determinant(), exponentInIntegral);

        // interpolatedDet += basisFunctions[3 * i + j] * pow( metrics[elemNodeTags[0][3 * iElm + j]].determinant(), exponentInIntegral);
        interpolatedDet += basisFunctions[3 * i + j] * metrics[elemNodeTags[0][3 * iElm + j]].determinant();
      }

      interpolatedDet = pow(interpolatedDet, exponentInIntegral);

      I += weights[i] * det[iElm * nQuad + i] * interpolatedDet;
      area += weights[i] * det[iElm * nQuad + i];
    }
  }

  feInfo("Computed integral   I = %1.5e", I);
  feInfo("Computed area = %1.5e", area);

  // const std::vector<Vertex> &meshVertices = _newRecovery->_mesh->getVertices();
  // for(auto &v : meshVertices) {
  //   MetricType M = metrics[&v];
  //   double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), exponentForDeterminant);
  //   M *= factor;
  //   metrics[&v] = M;
  // }

  for(size_t i = 0; i < nodeTags.size(); i++) {
    MetricType &M = metrics[nodeTags[i]];
    dim = 2.;
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), exponentForDeterminant);
    M *= factor;
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

template <class MetricType>
void feMetric::drawEllipsoids(const std::string &posFile, 
  const std::map<const Vertex*, MetricType> &metrics, double sizeFactor, int nPoints)
{
  FILE *f = fopen(posFile.data(), "w");
  fprintf(f, "View \"%s\"{\n", posFile.data());

  double factor = sizeFactor;
  int nt = nPoints;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  for(auto pair : metrics) {
    const Vertex *v = pair.first;
    MetricType m = pair.second;
    getEllipsePoints(factor * m(0, 0), factor * 2.0 * m(0, 1), factor * m(1, 1), v->x(), v->y(),
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

template <class MetricType>
void feMetric::drawEllipsoids(const std::string &posFile, 
  std::map<int, MetricType> &metrics, const std::vector<std::size_t> &nodeTags,
  const std::vector<double> &coord, double sizeFactor, int nPoints)
{
  FILE *f = fopen(posFile.data(), "w");
  fprintf(f, "View \"%s\"{\n", posFile.data());

  double factor = sizeFactor;
  int nt = nPoints;
  std::vector<double> xP(nt, 0.);
  std::vector<double> yP(nt, 0.);

  for(size_t i = 0; i < nodeTags.size(); i++) {
    MetricType &m = metrics[nodeTags[i]];
    const double x = coord[3*i+0];
    const double y = coord[3*i+1];
    getEllipsePoints(factor * m(0, 0), factor * 2.0 * m(0, 1), factor * m(1, 1), x, y, xP, yP);
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

feStatus feMetric::computeMetricsHechtKuate()
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

  return FE_STATUS_OK;
}

feStatus feMetric::computeMetricsLogSimplex()
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

  if(nodeTags.size() == 0) {
    feWarning("nodeTags is empty");
    exit(-1);
  }

  std::vector<Vertex> &meshVertices = _recovery->_mesh->getVertices();

  // Create a nodeTags to feVertex->tag map (brute-force for now)

  // On ne peut pas créer le map pour les maillages aniso intermédiaires, puisqu'il n'y a pas de
  // feMesh associé. Il faut recréer le maillage dans la boucle d'intégration pour pouvoir
  // interpoler les métriques. C'est-à-dire : il faut lancer avec -nLoop 0 -nAdapt N, puisque si
  // -nLoop est autre que 0 on boucle sur le maillage aniso avec la même reconstruction.
  double tol = 1e-5;
  _v2n.clear();
  _n2v.clear();
  for(auto &v : meshVertices) {
    _v2n[&v] = -1;
  }
  for(int i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    for(auto &v : meshVertices) {
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol) {
        _v2n[&v] = nodeTags[i];
        _n2v[nodeTags[i]] = &v;
      }
    }
  }
  // Check
  for(auto &v : meshVertices) {
    if(_v2n[&v] == -1) {
      feInfo("No Gmsh tag was associated to vertex %f - %f", v->x(), v->y());
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
  double lambdaMax = 1. / (lMin * lMin);
  double lambdaMin = 1. / (lMax * lMax);

  // Parameters for the log-simplex metric computation
  linearProblem myLP;
  Eigen::Matrix2d Q, V, D = Eigen::Matrix2d::Identity();
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;
  bool res;
  int maxIter = 50;
  tol = 1e-2; // Tolerance on the difference Q-Qprev
  int nThetaPerQuadrant = 250;
  int numIter;
  int size = 2 * 4 * nThetaPerQuadrant; // 2 coordonnees * 4 quadrants * nPhi angles/quadrant
  myLP.lvl1.resize(size, 0.);
  myLP.constraints.resize(size, 0.);

  // Set generic LP solver data
  myLP.problem.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
  myLP.problem.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
  DSVector dummycol(0);
  myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));
  myLP.problem.addColReal(LPCol(0.0, dummycol, infinity, -infinity));
  myLP.problem.addColReal(LPCol(1.0, dummycol, infinity, -infinity));

  myLP.row.setMax(3);
  myLP.prim.reDim(3);
  myLP.lprow = LPRow(3);
  myLP.lprow.setRhs(infinity);
  myLP.lprowset = LPRowSet(4 * nThetaPerQuadrant, 3);

  for(int i = 0; i < 4 * nThetaPerQuadrant; ++i) {
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

    res = computeMetricLogSimplexCurved(xpos, C, S, _recovery, Q, maxIter, nThetaPerQuadrant, tol,
                                        numIter, myLP);
    if(res) {
      feInfo("Metric (%+-1.5e - %+-1.5e - %+-1.5e) found in %2d iterations (log-simplex method) "
             "(vertex %d/%d) - Target vertices = %d",
             Q(0, 0), Q(1, 0), Q(1, 1), numIter, i, nodeTags.size(), _options.nTargetVertices);

      // Bound the eigenvalues of Q
      es.compute(Q, true);
      ev = es.eigenvalues();
      D(0, 0) = fmin(lambdaMax, fmax(lambdaMin, ev(0).real()));
      D(1, 1) = fmin(lambdaMax, fmax(lambdaMin, ev(1).real()));
      Q = es.eigenvectors().real() * D * es.eigenvectors().transpose().real();

    } else {
      feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e). Setting Q = lmax*I.", x, y);
      feInfo("C = %f - S = %f", C, S);
      Q = 1. / (lMax * lMax) * Eigen::Matrix2d::Identity();
      exit(-1);
    }

    SMetric3 M;
    M.set_m11(Q(0, 0));
    M.set_m21(Q(0, 1));
    M.set_m22(Q(1, 1));

    metrics[i] = M;
    _metricsOnGmshModel[nodeTags[i]] = M;

    // Plot the metrics raw metrics before scaling
    factor = 100.;
    getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[0], yP[0], 0., 1, 1);
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

    double dx = 2. / (NN - 1);
    double I = round((y + 1.) / dx) + 1;
    I = NN + 1 - I;
    double J = round((x + 1.) / dx) + 1;

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];
    fprintf(fMatlab, "%d %d %f %f %f %f\n", (int)I, (int)J, M(0, 0), M(0, 1), M(1, 0), M(1, 1));
  }
  fclose(fMatlab);

  if(_options.nTargetVertices > 0) {
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
        fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile5, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  if(_options.enableGradation) {
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
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];

    Eigen::Matrix2d m;
    m(0, 0) = M(0, 0);
    m(0, 1) = M(0, 1);
    m(1, 0) = M(1, 0);
    m(1, 1) = M(1, 1);
    // ////////////////////////////////////////////
    // m(0,0) = 10.456 + 2.123 * x*x + 4.6789 * y;
    // m(0,1) = 0.;
    // m(1,0) = 0.;
    // m(1,1) = 10.456 + 2.123 * x*x + 4.6789 * y;
    // ////////////////////////////////////////////
    _metricsOnGmshModel_eigen[nodeTags[i]] = m;

    std::vector<double> vMetric(9);
    std::vector<double> vRecovery(1);

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
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

  gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           recoveryData);

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);
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

  if(false) {
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

    for(int ielm = 0; ielm < _recovery->_mesh->_elements.size(); ++ielm) {
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
      double x12[2] = {(p0[0] + p1[0]) / 2. + a[0], (p0[1] + p1[1]) / 2. + a[1]};
      // interpolateMetricP1(x12, logM11, M0, M1, M2, sumduda1M, sumduda2M, false);

      // for(int it = 0; it <= 10; ++it){
      //   double tt = (double) it / 10.;
      //   double pos[2] = {p0[0] + gamma[0]*tt + 4.*tt*(1.-tt)*a[0], p0[1] + gamma[1]*tt
      //   + 4.*tt*(1.-tt)*a[1]};

      //   interpolateMetricAndDerivativeOnP2Edge(tt, logM11, M0, M1, sumduda1M, sumduda2M, myM,
      //   dMda1, dMda2);

      //   factor = 100.;
      //   bool res = getEllipsePoints(factor * myM(0, 0), factor * 2.0 * myM(0, 1), factor * myM(1,
      //   1), pos[0], pos[1], xP, yP); if(res){
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
      //   res = getEllipsePoints(factor * dMda1(0, 0), factor * 2.0 * dMda1(0, 1), factor *
      //   dMda1(1, 1), pos[0], pos[1], xP, yP); if(res){
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
      //   res = getEllipsePoints(factor * dMda2(0, 0), factor * 2.0 * dMda2(0, 1), factor *
      //   dMda2(1, 1), pos[0], pos[1], xP, yP); if(res){
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
      for(int i = 0; i <= ni; ++i) {
        for(int j = 0; j <= nj; ++j) {
          double rr = (double)i / (double)ni;
          double ss = (double)j / (double)nj;
          if(rr + ss <= 1.) {
            double pos[2] = {(1 - rr - ss) * v1->x() + rr * v2->x() + ss * v3->x(),
                             (1 - rr - ss) * v1->y() + rr * v2->y() + ss * v3->y()};

            interpolateMetricP1(pos, myM, sumduda1M, sumduda2M);

            // myM = myM.exp();
            factor = 100.;
            bool res = getEllipsePoints(factor * myM(0, 0), factor * 2.0 * myM(0, 1),
                                        factor * myM(1, 1), pos[0], pos[1], xP, yP);
            if(res) {
              for(int j = 0; j < nt; ++j) {
                if(j != nt - 1) {
                  fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                          xP[j + 1], yP[j + 1], 0., 1, 1);
                } else {
                  fprintf(f, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                          xP[0], yP[0], 0., 1, 1);
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

  return FE_STATUS_OK;

#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to create a metric field and SoPlex is required "
         "to compute the metric tensors.\n");
#endif
}

// The metric field is stored in the view _metricViewTag in the
// gmsh model given in metricOptions.gmshModel
feStatus feMetric::computeMetricsExtremeSizesOnly()
{
#if defined(HAVE_GMSH)
  if(!_options.isGmshModelReady) {
    return feErrorMsg(FE_STATUS_ERROR, "No Gmsh model available to compute a metric field.\n"
           "Create a Gmsh geometry and mesh it first.\n");
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

  if(nodeTags.size() == 0) {
    return feErrorMsg(FE_STATUS_ERROR, "nodeTags is empty");
  }

  std::vector<Vertex> &meshVertices = _recovery->_mesh->getVertices();

  // Create a nodeTags to feVertex->tag map (brute-force for now)

  // On ne peut pas créer le map pour les maillages aniso intermédiaires, puisqu'il n'y a pas de
  // feMesh associé. Il faut recréer le maillage dans la boucle d'intégration pour pouvoir
  // interpoler les métriques. C'est-à-dire : il faut lancer avec -nLoop 0 -nAdapt N, puisque si
  // -nLoop est autre que 0 on boucle sur le maillage aniso avec la même reconstruction.
  double tol = 1e-5;
  _v2n.clear();
  _n2v.clear();
  for(auto &v : meshVertices) {
    _v2n[&v] = -1;
  }
  for(int i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    for(auto &v : meshVertices) {
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol) {
        _v2n[&v] = nodeTags[i];
        _n2v[nodeTags[i]] = &v;
      }
    }
  }
  // Check
  for(auto &v : meshVertices) {
    if(_v2n[&v] == -1) {
      return feErrorMsg(FE_STATUS_ERROR, "No Gmsh tag was associated to vertex %f - %f", v.x(), v.y());
    }
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

  FILE *graBeforeSmoothing = fopen("graBeforeSmoothing.pos", "w");
  FILE *isoBeforeSmoothing = fopen("isoBeforeSmoothing.pos", "w");
  fprintf(graBeforeSmoothing, "View \" graBeforeSmoothing \"{\n");
  fprintf(isoBeforeSmoothing, "View \" isoBeforeSmoothing \"{\n");

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    xpos[0] = x;
    xpos[1] = y;

    switch(_options.directionFieldFromDerivativesOfOrder) {
      case 1:
        computeDirectionFieldFromGradient(xpos, C, S, 1e-8, _recovery, graBeforeSmoothing,
                                          isoBeforeSmoothing);
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

  fprintf(graBeforeSmoothing, "};");
  fclose(graBeforeSmoothing);
  fprintf(isoBeforeSmoothing, "};");
  fclose(isoBeforeSmoothing);

  // int smoothMaxIter = 100;
  // double smoothTol = 0.8;
  // smoothDirections(COS, SIN, smoothMaxIter, smoothTol);

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
  double lambdaMax = 1. / (lMin * lMin);
  double lambdaMin = 1. / (lMax * lMax);

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];

    xpos[0] = x;
    xpos[1] = y;

    double C = COS[nodeTags[i]];
    double S = SIN[nodeTags[i]];

    double hIso, hGrad;
    switch(deg) {
      case 2: {
        const double derivativeAlongGradient = fabs(dttt(xpos, C, S, _recovery, 0));
        const double derivativeAlongIsoline = fabs(dttt(xpos, C, S, _recovery, 1));

        hIso = pow(6.0 * eps / derivativeAlongIsoline, 0.3333);
        hGrad = pow(6.0 * eps / derivativeAlongGradient, 0.3333);

        break;
      }
      default: {
        return feErrorMsg(FE_STATUS_ERROR, "In feMetric : Error - No metric computation scheme for deg != 2\n");
      }
    }

    Eigen::Matrix2d Q, V, D = Eigen::Matrix2d::Identity();
    Eigen::EigenSolver<Eigen::Matrix2d> es;
    Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;

    double lambdaIso = 1. / (hIso * hIso);
    double lambdaGrad = 1. / (hGrad * hGrad);

    V(0, 0) = C;
    V(1, 0) = S;

    V(0, 1) = -S;
    V(1, 1) = C;

    D(0, 0) = fmin(lambdaMax, fmax(lambdaMin, lambdaGrad));
    D(1, 1) = fmin(lambdaMax, fmax(lambdaMin, lambdaIso));

    Q = V * D * V.transpose();

    SMetric3 M;
    M.set_m11(Q(0, 0));
    M.set_m21(Q(0, 1));
    M.set_m22(Q(1, 1));

    metrics[i] = M;
    _metricsOnGmshModel[nodeTags[i]] = M;

    // Plot the "raw" metrics before scaling
    factor = 100.;
    getEllipsePoints(factor * M(0, 0), factor * 2.0 * M(0, 1), factor * M(1, 1), x, y, xP, yP);
    for(int j = 0; j < nt; ++j) {
      if(j != nt - 1) {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile3, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  if(_options.enableGradation) {
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
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[j + 1], yP[j + 1], 0., 1, 1);
      } else {
        fprintf(debugFile6, "SL(%.16g,%.16g,%.16g,%.16g,%.16g,%.16g){%u, %u};\n", xP[j], yP[j], 0.,
                xP[0], yP[0], 0., 1, 1);
      }
    }
  }

  for(size_t i = 0; i < nodeTags.size(); i++) {

    SMetric3 M = _metricsOnGmshModel[nodeTags[i]];

    Eigen::Matrix2d m;
    m(0, 0) = M(0, 0);
    m(0, 1) = M(0, 1);
    m(1, 0) = M(1, 0);
    m(1, 1) = M(1, 1);
    // ////////////////////////////////////////////
    // m(0,0) = 10.456 + 2.123 * x + 4.6789 * y;
    // m(0,1) = 0.;
    // m(1,0) = 0.;
    // m(1,1) = 10.456 + 2.123 * x + 4.6789 * y;
    // ////////////////////////////////////////////
    _metricsOnGmshModel_eigen[nodeTags[i]] = m;

    std::vector<double> vMetric(9);
    std::vector<double> vRecovery(1);

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
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

  gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           recoveryData);

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);
  gmsh::view::write(recoveryViewTag, _options.recoveryName);

  return FE_STATUS_OK;

#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to create a metric field and SoPlex is required "
         "to compute the metric tensors.\n");
#endif
}

template<class MetricType>
static SMetric3 convert2metric3(const MetricType &other)
{
  SMetric3 M(1.0);
  M(0,0) = other(0,0);
  M(1,0) = other(1,0);
  M(0,1) = other(0,1);
  M(1,1) = other(1,1);
  return M;
}

//
// Metric tensors for P1 adaptation (Alauzet & Loseille)
//
feStatus feMetric::computeMetricsP1()
{
#if defined(HAVE_GMSH)
  if(!_options.isGmshModelReady) {
    return feErrorMsg(FE_STATUS_ERROR, "No Gmsh model available to compute a metric field.\n"
           "Create a Gmsh geometry and mesh it first.\n");
  }

  feInfo("Computing metrics on Gmsh model %s", _options.modelForMetric.c_str());
  gmsh::model::setCurrent(_options.modelForMetric);

  // Get the nodes from the Gmsh model
  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  std::vector<double> parametricCoord;
  int dim = -1;
  int tag = -1;
  int includeBoundary = false;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, includeBoundary, false);
  // getNodesForPhysicalGroup( ... )

  if(nodeTags.size() == 0) {
    return feErrorMsg(FE_STATUS_ERROR, "Trying to compute metric tensors on a mesh with 0 nodes");
  }

  // Min and max eigenvalues based on sizes
  double lambdaMax = 1. / (_options.hMin * _options.hMin);
  double lambdaMin = 1. / (_options.hMax * _options.hMax);

  // Create a view which will contain the metric field
  _metricViewTag = gmsh::view::add(":metric", 17);
  feInfo("_metricViewTag = %d", _metricViewTag);
  int recoveryViewTag = gmsh::view::add("recovery", 21);
  std::vector<std::vector<double> > metricData;
  std::vector<std::vector<double> > recoveryData;

  double fxx, fxy, fyx, fyy;
  MetricTensor H2;
  Eigen::Matrix2d H, D = Eigen::Matrix2d::Identity();
  Eigen::EigenSolver<Eigen::Matrix2d> es;
  Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType ev;

  double x[2];
  // Compute hessian at vertices
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    fxx = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 2, x);
    fxy = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 3, x);
    fyx = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 4, x);
    fyy = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 5, x);

    H(0,0) = fxx;
    H(0,1) = (fxy+fyx)/2.;
    H(1,0) = (fxy+fyx)/2.;
    H(1,1) = fyy;

    // Bound the eigenvalues of H
    es.compute(H, true);
    ev = es.eigenvalues();
    D(0, 0) = fmin(lambdaMax, fmax(lambdaMin, fabs(ev(0).real()) ));
    D(1, 1) = fmin(lambdaMax, fmax(lambdaMin, fabs(ev(1).real()) ));
    H = es.eigenvectors().real() * D * es.eigenvectors().transpose().real();

    _metricsP1AtNodeTags[nodeTags[i]] = H;

    H2(0,0) = fxx;
    H2(0,1) = (fxy+fyx)/2.;
    H2(1,0) = (fxy+fyx)/2.;
    H2(1,1) = fyy;
    _metricTensorAtNodetags[nodeTags[i]] = H2.boundEigenvaluesOfAbs(lambdaMin, lambdaMax);
  }

  // drawEllipsoids("ellipses.pos", _metricsP1AtNodeTags, nodeTags, coord, 100, 30);
  drawEllipsoids("ellipses.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);

  double N = (double) _options.nTargetVertices;
  double p = _options.LpNorm;
  feInfo("Targetting %f vertices in L%f norm", N, p);
  double n = 2.; // Space dimension
  double exponentInIntegral =       p / (2.*p + n);
  double exponentForDeterminant = -1. / (2.*p + n);

  // metricScalingFromGmshSubstitute(_metricsP1AtNodeTags, nodeTags, exponentInIntegral, exponentForDeterminant);
  metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, exponentInIntegral, exponentForDeterminant);

  // drawEllipsoids("afterScaling.pos", _metricsP1AtNodeTags, nodeTags, coord, 100, 30);
  drawEllipsoids("afterScaling.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);

  // Transfer from one map to the other to compute gradation (FIXME)
  for(size_t i = 0; i < nodeTags.size(); i++) {
    // _metricsP1[nodeTags[i]] = convert2metric3(_metricsP1AtNodeTags[nodeTags[i]]);
    _metricsP1[nodeTags[i]] = convert2metric3(_metricTensorAtNodetags[nodeTags[i]]);
  }

  if(_options.enableGradation) {
    gradationMetriques(_options.gradation, 200, coord, _metricsP1);
  }

  // Transfer back
  for(size_t i = 0; i < nodeTags.size(); i++) {
    // _metricsP1AtNodeTags[nodeTags[i]](0,0) = _metricsP1[nodeTags[i]](0,0);
    // _metricsP1AtNodeTags[nodeTags[i]](1,0) = _metricsP1[nodeTags[i]](1,0);
    // _metricsP1AtNodeTags[nodeTags[i]](0,1) = _metricsP1[nodeTags[i]](0,1);
    // _metricsP1AtNodeTags[nodeTags[i]](1,1) = _metricsP1[nodeTags[i]](1,1);
    _metricTensorAtNodetags[nodeTags[i]](0,0) = _metricsP1[nodeTags[i]](0,0);
    _metricTensorAtNodetags[nodeTags[i]](1,0) = _metricsP1[nodeTags[i]](1,0);
    _metricTensorAtNodetags[nodeTags[i]](0,1) = _metricsP1[nodeTags[i]](0,1);
    _metricTensorAtNodetags[nodeTags[i]](1,1) = _metricsP1[nodeTags[i]](1,1);
  }

  // drawEllipsoids("afterGradation.pos", _metricsP1AtNodeTags, nodeTags, coord, 100, 30);
  drawEllipsoids("afterGradation.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);

  std::vector<double> vMetric(9);
  std::vector<double> vRecovery(1);
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    // For interpolation in curved mesh
    // _metricsOnGmshModel_eigen[nodeTags[i]] = m;

    // Eigen::Matrix2d M = _metricsP1AtNodeTags[nodeTags[i]];
    MetricTensor &M = _metricTensorAtNodetags[nodeTags[i]];

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.0;
    vRecovery[0] = _newRecovery->evaluateRecovery(PPR::RECOVERY, 0, x);;

    metricData.push_back(vMetric);
    recoveryData.push_back(vRecovery);
  }

  gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           recoveryData);

  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);
  gmsh::view::write(recoveryViewTag, _options.recoveryName);

  return FE_STATUS_OK;

#else
  return feErrorMsg(FE_STATUS_ERROR, "Gmsh is required to create a metric field and SoPlex is required "
         "to compute the metric tensors.\n");
#endif
}

feStatus feMetric::computeMetrics()
{
  switch(_options.method) {
    case adaptationMethod::ANISO_PN:
      return computeMetricsHechtKuate();

    case adaptationMethod::ANISO_P1:
      return computeMetricsP1();

    case adaptationMethod::CURVED_LS:
      return computeMetricsLogSimplex();

    case adaptationMethod::CURVED_EXTREME_SIZES:
      return computeMetricsExtremeSizesOnly();
  }
  return FE_STATUS_OK;
}