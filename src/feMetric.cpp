
#include "ellipseToolbox.h"
#include "fullMatrix.h"
#include "feMessage.h"
#include "feMetric.h"
#include "feMetricTools.h"
#include "../contrib/root_finder/root_finder.hpp"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#if defined(HAVE_SOPLEX)
#include "soplex.h"
using namespace soplex;
#endif

extern int FE_VERBOSE;

// #define ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE
// #define TEST_ANALYTIC_METRIC

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions)
  : _options(metricOptions), _recovery(recovery)
{
  _options.polynomialDegree = _recovery->getDegreeSolution();
}

feMetric::feMetric(std::vector<feNewRecovery *> &recoveredFields, feMetricOptions metricOptions)
  : _options(metricOptions), _recoveredFields(recoveredFields)
{
  // Min and max eigenvalues based on sizes
  _lambdaMax = 1. / (_options.hMin * _options.hMin);
  _lambdaMin = 1. / (_options.hMax * _options.hMax);

  _options.polynomialDegree = _recoveredFields[0]->getDegreeSolution();
}

feStatus feMetric::createVertex2NodeMap(std::vector<std::size_t> &nodeTags,
                                        std::vector<double> &coord)
{
#if defined(HAVE_GMSH)

  // Get the nodeTags and coordinates from the Gmsh model
  std::vector<double> parametricCoord;
  int dim = -1;
  int tag = -1;
  int includeBoundary = false;
  gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, dim, tag, includeBoundary, false);

  // Create a nodeTags to feVertex->tag map (brute-force for now)
  feMesh *mesh;
  if(_recoveredFields[0] != nullptr)
    mesh = _recoveredFields[0]->_mesh;
  else if(_recovery != nullptr)
    mesh = _recovery->_mesh;
  else
    return feErrorMsg(FE_STATUS_ERROR, "No valid recovery");

  std::vector<Vertex> &meshVertices = mesh->getVertices();

  double tol = 1e-12;
  _v2n.clear();
  _n2v.clear();

  _sequentialTag2nodeTagVec.resize(meshVertices.size());

  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    for(auto &v : meshVertices) {
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol) {
        _v2n[&v] = nodeTags[i];
        _n2v[nodeTags[i]] = &v;
        _nodeTag2sequentialTag[nodeTags[i]] = mesh->_verticesMap[v.getTag()];
        _sequentialTag2nodeTag[mesh->_verticesMap[v.getTag()]] = nodeTags[i];
        _sequentialTag2nodeTagVec[mesh->_verticesMap[v.getTag()]] = nodeTags[i];
        break;
      }
    }
  }
  // Debug check (not valid when computing metrics on P1 mesh from derivatives on a P2 mesh)
  // for(auto &v : meshVertices) {
  //   auto it = _v2n.find(&v);
  //   if(it == _v2n.end())
  //     return feErrorMsg(FE_STATUS_ERROR, "No Gmsh tag was associated to vertex %f - %f", v.x(),
  //                       v.y());
  // }
#else
  UNUSED(nodeTags, coord);
#endif
  return FE_STATUS_OK;
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is
// performed on the Gmsh subsitute (gmshModel). Using gmsh api because it's easier.
// FIXME: parallelize
template <class MetricType>
void feMetric::metricScalingFromGmshSubstitute(std::map<int, MetricType> &metrics,
                                               const std::vector<size_t> &nodeTags,
                                               const std::vector<double> &/*coord*/,
                                               double exponentInIntegral,
                                               double exponentForDeterminant)
{
#if defined(HAVE_GMSH)

  int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  if(elementTypes.size() > 1) {
    feErrorMsg(FE_STATUS_ERROR, "Gmsh mesh has more than one element type."
                                "We expect only P1 or P2 triangles in 2D entities");
    exit(-1);
  }

  int elementType = elementTypes[0];

  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  if(_nVerticesPerElmOnBackmesh == 3 && nVerticesPerElement == 6) {
    // Linear backmesh with a P2 feMesh : get P1 basis
    int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
    gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents,
                                         basisFunctions, numOrientations);
  }

  // Get the jacobians
  std::vector<double> jac, detJac, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, detJac, pts);

  // Compute integral of det^exponent
  double I = 0.0;
  size_t nQuad = weights.size();
  size_t nElm = elementTags[0].size();

  double N = (double)_options.nTargetVertices;
  double dim = 2.;

#define MAX_VERTICES_BACKMESH 6

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    MetricType M_interpolated;
    double xsi[2];
    int tags[MAX_VERTICES_BACKMESH];

#if defined(HAVE_OMP)
#pragma omp for reduction(+ : I)
#endif
    for(size_t iElm = 0; iElm < nElm; iElm++) {

      for(int i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
        tags[i] = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + i];
      }

      MetricTensor logM0 = metrics.at(tags[0]).log();
      MetricTensor logM1 = metrics.at(tags[1]).log();
      MetricTensor logM2 = metrics.at(tags[2]).log();

      for(size_t i = 0; i < nQuad; i++) {
        xsi[0] = localCoord[3 * i + 0];
        xsi[1] = localCoord[3 * i + 1];

        // Coordonnees physiques
        // double xphys = pts[iElm * 3 * nQuad + 3 * i + 0];
        // double yphys = pts[iElm * 3 * nQuad + 3 * i + 1];

        if(_nVerticesPerElmOnBackmesh == 3) {
          logEuclidianP1Interpolation(xsi, logM0, logM1, logM2, M_interpolated);
        }

        if(_nVerticesPerElmOnBackmesh == 6) {
          MetricTensor logM3 = metrics.at(tags[3]).log();
          MetricTensor logM4 = metrics.at(tags[4]).log();
          MetricTensor logM5 = metrics.at(tags[5]).log();
          logEuclidianP2Interpolation(xsi, logM0, logM1, logM2, logM3, logM4, logM5,
                                      M_interpolated);
        }

        double interpolatedDet = M_interpolated.determinant();

        if(interpolatedDet <= 0.) {
          feInfo("Negative metric determinant:");
          M_interpolated.print();
          exit(-1);
        }

        I += weights[i] * detJac[iElm * nQuad + i] * pow(interpolatedDet, exponentInIntegral);
      }
    }
  }

  feInfoCond(FE_VERBOSE > 0, "Scaling : Integral of determinant = %+-1.16e", I);

  // Apply scaling
  for(size_t i = 0; i < nodeTags.size(); i++) {
    MetricType &M = metrics[nodeTags[i]];
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), exponentForDeterminant);
    M *= factor;
  }
#else
  UNUSED(metrics, nodeTags, exponentInIntegral, exponentForDeterminant);
#endif
}

//////////////////////////////////////////////////////////////
// Begining of metric interpolation convergence tests
//
double matNorm2(const MetricTensor &m1, const MetricTensor &m2)
{
  double sqr = 0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      sqr += (m1(i, j) - m2(i, j)) * (m1(i, j) - m2(i, j));
    }
  }
  return sqrt(sqr);
}

inline double matfun(double x, double y) { return 1. + x * x * x + y * y; }

MetricTensor analyticMetric(double x, double y)
{
  MetricTensor res(1.0);
  SVector3 EX(1., 0., 0.);
  SVector3 EY(0., 1., 0.);

  SVector3 e1 = EX * cos(M_PI * x) + EY * sin(M_PI * x);
  SVector3 e2 = EX * (-sin(M_PI * x)) + EY * cos(M_PI * x);
  SVector3 e3(0., 0., 1.);
  double h1 = 1 + x * y;
  double h2 = 5 + x * y * y;
  SMetric3 tmp(1. / (h1 * h1), 1. / (h2 * h2), 1., e1, e2, e3);
  res(0, 0) = tmp(0, 0);
  res(0, 1) = tmp(0, 1);
  res(1, 0) = tmp(1, 0);
  res(1, 1) = tmp(1, 1);
  return res;
}

void feMetric::interpolationTest(const std::vector<size_t> &nodeTags, std::vector<double> &coord)
{
#if defined(HAVE_GMSH)

  int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int elementOrder = getGeometricInterpolantDegree(_recoveredFields[0]->_cnc->getInterpolant());
  feInfo("Mesh elements order: %d", elementOrder);
  int elementType = gmsh::model::mesh::getElementType("Triangle", elementOrder);

  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  // Get the mesh elements
  // std::vector<int> elementTypes;
  // std::vector<std::vector<std::size_t> > elementTags;
  // std::vector<std::vector<std::size_t> > elemNodeTags;
  // gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);
  std::vector<std::size_t> elementTags;
  std::vector<std::size_t> elemNodeTags;
  gmsh::model::mesh::getElementsByType(elementType, elementTags, elemNodeTags);

  if(elementTags.empty()) {
    feErrorMsg(FE_STATUS_ERROR, "There are no P2 elements in the Gmsh model. It's probably because "
                                "the Gmsh model is still P1 for MMG.");
    exit(-1);
  }
  feInfo("Number of elements in the Gmsh model = %d", elementTags.size());

  // Get the jacobians
  std::vector<double> jac, det, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, det, pts);

  // Compute integral of det^exponent
  double area = 0.0;
  int nQuad = weights.size();

  feInfo("Size of localCoord = %d", localCoord.size());
  feInfo("nQuad = %d", nQuad);

  double error = 0.;

  double xsi[2];

  double scale = 10000;

  FILE *myfile = fopen("checkCoord.pos", "w");
  FILE *myfile2 = fopen("metricExact.pos", "w");
  FILE *myfile3 = fopen("metricInterp.pos", "w");
  FILE *myfile4 = fopen("metricSommets.pos", "w");
  fprintf(myfile, "View\"coord\"{\n");
  fprintf(myfile2, "View\"metricExact\"{\n");
  fprintf(myfile3, "View\"metricInterp\"{\n");
  fprintf(myfile4, "View\"metricSommets\"{\n");

  for(size_t iElm = 0; iElm < elementTags.size(); iElm++) {
    int v0 = elemNodeTags[nVerticesPerElement * iElm + 0] - 1;
    int v1 = elemNodeTags[nVerticesPerElement * iElm + 1] - 1;
    int v2 = elemNodeTags[nVerticesPerElement * iElm + 2] - 1;
    int v3 = elemNodeTags[nVerticesPerElement * iElm + 3] - 1;
    int v4 = elemNodeTags[nVerticesPerElement * iElm + 4] - 1;
    int v5 = elemNodeTags[nVerticesPerElement * iElm + 5] - 1;

    double x0 = coord[3 * v0];
    double y0 = coord[3 * v0 + 1];
    double x1 = coord[3 * v1];
    double y1 = coord[3 * v1 + 1];
    double x2 = coord[3 * v2];
    double y2 = coord[3 * v2 + 1];
    double x3 = coord[3 * v3];
    double y3 = coord[3 * v3 + 1];
    double x4 = coord[3 * v4];
    double y4 = coord[3 * v4 + 1];
    double x5 = coord[3 * v5];
    double y5 = coord[3 * v5 + 1];

    if(iElm == 0) {
      fprintf(myfile, "SP(%f,%f,0.){0.};\n", x0, y0);
      fprintf(myfile, "SP(%f,%f,0.){1.};\n", x1, y1);
      fprintf(myfile, "SP(%f,%f,0.){2.};\n", x2, y2);
      fprintf(myfile, "SP(%f,%f,0.){3.};\n", x3, y3);
      fprintf(myfile, "SP(%f,%f,0.){4.};\n", x4, y4);
      fprintf(myfile, "SP(%f,%f,0.){5.};\n", x5, y5);
    }

    double X[6] = {x0, x1, x2, x3, x4, x5};
    double Y[6] = {y0, y1, y2, y3, y4, y5};
    double U[6] = {matfun(x0, y0), matfun(x1, y1), matfun(x2, y2),
                   matfun(x3, y3), matfun(x4, y4), matfun(x5, y5)};

    // MetricTensor M0(1.0); M0(0,0) = matfun(x0,y0); M0(1,1) = matfun(x0,y0);
    // MetricTensor M1(1.0); M1(0,0) = matfun(x1,y1); M1(1,1) = matfun(x1,y1);
    // MetricTensor M2(1.0); M2(0,0) = matfun(x2,y2); M2(1,1) = matfun(x2,y2);
    // MetricTensor M3(1.0); M3(0,0) = matfun(x3,y3); M3(1,1) = matfun(x3,y3);
    // MetricTensor M4(1.0); M4(0,0) = matfun(x4,y4); M4(1,1) = matfun(x4,y4);
    // MetricTensor M5(1.0); M5(0,0) = matfun(x5,y5); M5(1,1) = matfun(x5,y5);

    MetricTensor M0 = analyticMetric(x0, y0);
    MetricTensor M1 = analyticMetric(x1, y1);
    MetricTensor M2 = analyticMetric(x2, y2);
    MetricTensor M3 = analyticMetric(x3, y3);
    MetricTensor M4 = analyticMetric(x4, y4);
    MetricTensor M5 = analyticMetric(x5, y5);

    if(iElm == 0) {
      double POS0[2] = {x0, y0};
      double POS1[2] = {x1, y1};
      double POS2[2] = {x2, y2};
      double POS3[2] = {x3, y3};
      double POS4[2] = {x4, y4};
      double POS5[2] = {x5, y5};
      drawSingleEllipse(myfile4, POS0, M0, scale, 100);
      drawSingleEllipse(myfile4, POS1, M1, scale, 100);
      drawSingleEllipse(myfile4, POS2, M2, scale, 100);
      drawSingleEllipse(myfile4, POS3, M3, scale, 100);
      drawSingleEllipse(myfile4, POS4, M4, scale, 100);
      drawSingleEllipse(myfile4, POS5, M5, scale, 100);
    }

    for(int i = 0; i < nQuad; i++) {
      // Interpolate metric at quad node
      MetricTensor Mk(1.0);
      xsi[0] = localCoord[3 * i + 0];
      xsi[1] = localCoord[3 * i + 1];
      // classicalP1Interpolation(xsi, M0, M1, M2, Mk);
      // logEuclidianP1Interpolation(xsi, M0, M1, M2, Mk);
      logEuclidianP2Interpolation(xsi, M0, M1, M2, M3, M4, M5, Mk);

      // Exact metric at quad node
      double xphys = 0.;
      double yphys = 0.;
      double uh = 0.;
      double sum = 0.;
      for(int ii = 0; ii < 6; ++ii) {
        sum += basisFunctions[nVerticesPerElement * i + ii];
        xphys += basisFunctions[nVerticesPerElement * i + ii] * X[ii];
        yphys += basisFunctions[nVerticesPerElement * i + ii] * Y[ii];
        uh += basisFunctions[nVerticesPerElement * i + ii] * U[ii];
      }

      // MetricTensor Mkexact(1.0);
      // Mkexact(0,0) = matfun(xphys,yphys);
      // Mkexact(1,1) = matfun(xphys,yphys);
      MetricTensor Mkexact = analyticMetric(xphys, yphys);

      if(iElm == 0) {
        double POS_QUAD[2] = {xphys, yphys};
        drawSingleEllipse(myfile2, POS_QUAD, Mkexact, scale, 100);
        drawSingleEllipse(myfile3, POS_QUAD, Mk, scale, 100);
      }

      // double uexact = matfun(xphys, yphys);

      double frobError = matNorm2(Mk, Mkexact);

      // Erreur sur l'interpolation des metriques
      error += frobError * frobError * weights[i] * det[iElm * nQuad + i];

      // Erreur (test de verification) sur l'interpolation d'un champ
      // error += (uh-uexact)*(uh-uexact) * weights[i] * det[iElm * nQuad + i];

      area += weights[i] * det[iElm * nQuad + i];
    }
  }

  fprintf(myfile, "};\n");
  fclose(myfile);
  fprintf(myfile2, "};\n");
  fclose(myfile2);
  fprintf(myfile3, "};\n");
  fclose(myfile3);
  fprintf(myfile4, "};\n");
  fclose(myfile4);

  error = sqrt(error);
  feInfo("Error is %f", error);

  // feInfo("Computed integral   I = %1.5e", I);
  feInfo("Computed area = %1.5e", area);

  FILE *file = fopen("checkInterpolation.pos", "w");
  fprintf(file, "View\"interpolation\"{\n");

  double x[2];

  // Apply scaling
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    // MetricTensor M(1.0);
    // M(0,0) = 1.0 + x[0]*x[0];
    // M(0,1) = 0.;
    // M(1,0) = 0.;
    // M(1,1) = 1.0 + x[0]*x[0];

    MetricTensor M = analyticMetric(x[0], x[1]);

    drawSingleEllipse(file, x, M, 500., 30);
  }

  fprintf(file, "};\n");
  fclose(file);
  _options.userValue = error;
#else
  UNUSED(nodeTags, coord);
#endif
}
//
// End of metric interpolation tests
///////////////////////////////////////////////////////////////////////

// Return a SMetric3 (Gmsh type for symmetric matrix) from a MetricType
// MetricType can be either an Eigen::Matrix2d or a MetricTensor (our custom type)
template <class MetricType> SMetric3 convert2metric3(const MetricType &other)
{
  SMetric3 M(1.0);
  M(0, 0) = other(0, 0);
  M(1, 0) = other(1, 0);
  M(0, 1) = other(0, 1);
  M(1, 1) = other(1, 1);
  return M;
}

Eigen::Matrix2d convert2eigenMatrix(const MetricTensor &other)
{
  Eigen::Matrix2d M;
  M(0, 0) = other(0, 0);
  M(1, 0) = other(1, 0);
  M(0, 1) = other(0, 1);
  M(1, 1) = other(1, 1);
  return M;
}

// Apply Alauzet's metric gradation
// FIXME: gradation is computed for SMetric3 metric for now, so
// need to transfer MetricTensors to SMetric3, then transfer back
void feMetric::applyGradation(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
  // Transfer from one map to the other to compute gradation
  for(size_t i = 0; i < nodeTags.size(); i++) {
    _smetric3AtNodetags[nodeTags[i]] = convert2metric3(_metricTensorAtNodetags[nodeTags[i]]);
  }

  gradationMetriques(_options.gradation, 200, coord, _smetric3AtNodetags);

  // Transfer back
  for(size_t i = 0; i < nodeTags.size(); i++) {
    _metricTensorAtNodetags[nodeTags[i]](0, 0) = _smetric3AtNodetags[nodeTags[i]](0, 0);
    _metricTensorAtNodetags[nodeTags[i]](1, 0) = _smetric3AtNodetags[nodeTags[i]](1, 0);
    _metricTensorAtNodetags[nodeTags[i]](0, 1) = _smetric3AtNodetags[nodeTags[i]](0, 1);
    _metricTensorAtNodetags[nodeTags[i]](1, 1) = _smetric3AtNodetags[nodeTags[i]](1, 1);
  }
}

void feMetric::writeMetricField(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
#if defined(HAVE_GMSH)

  std::vector<std::size_t> P1ToP2;
  size_t numBGMVertices = nodeTags.size();
  std::vector<std::size_t> nodeTagsP1;
  std::string initialMesh;
  gmsh::model::getCurrent(initialMesh);

  // If the background mesh elements are P2, create a P1 copy of the mesh
  // and export the metrics on the P1 mesh for MMG.
  if(_backmeshOrder > 1) {
    // Write a P1 .msh file to give MMG for aniso straight adaptation
    // For now copy the mesh and read it again
    std::string cmd = "cp " + _options.backgroundMeshfile + " tmp_P1_bgm.msh";
    int retvalue = system(cmd.data());
    if(retvalue != 0) {
      feErrorMsg(FE_STATUS_WRITE_ERROR, "Could not copy background mesh file! (return value = %d)",
        retvalue);
    }
    gmsh::open("tmp_P1_bgm.msh");

    gmsh::model::mesh::setOrder(1);
    gmsh::model::mesh::renumberNodes();

    // Brute-force identification of P1 vertices in P2 mesh
    double tol = 1e-8;
    std::vector<double> coordP1;
    std::vector<double> parametricCoord;
    int dim = -1;
    int tag = -1;
    int includeBoundary = false;
    gmsh::model::mesh::getNodes(nodeTagsP1, coordP1, parametricCoord, dim, tag, includeBoundary,
                                false);

    // P1 to P2 nodeTags map
    numBGMVertices = nodeTagsP1.size();
    P1ToP2.resize(numBGMVertices);

    for(size_t i = 0; i < numBGMVertices; ++i) {
      for(size_t j = 0; j < nodeTags.size(); ++j) {
        double xp1 = coordP1[3 * i + 0];
        double yp1 = coordP1[3 * i + 1];
        double xp2 = coord[3 * j + 0];
        double yp2 = coord[3 * j + 1];
        if(sqrt((xp1 - xp2) * (xp1 - xp2) + (yp1 - yp2) * (yp1 - yp2)) < tol) {
          P1ToP2[i] = nodeTags[j];
          break;
        }
      }
    }
  }

  std::string bgmToWrite;
  gmsh::model::getCurrent(bgmToWrite);

  // Create a view that contains the metric field
  _metricViewTag = gmsh::view::add(":metric");
  // int recoveryViewTag = gmsh::view::add(":recovery");
  std::vector<std::vector<double> > metricData;
  // std::vector<std::vector<double>> recoveryData;

  // double x[2];
  int tag;
  std::vector<double> vMetric(9);
  // std::vector<double> vRecovery(1);

  for(size_t i = 0; i < numBGMVertices; i++) {
    // x[0] = coord[3 * i + 0];
    // x[1] = coord[3 * i + 1];

    // Get metric at vertex (if P1 background mesh) or at corresponding P2 vertex
    tag = (_backmeshOrder > 1) ? P1ToP2[i] : nodeTags[i];

    MetricTensor &M = _metricTensorAtNodetags[tag];

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.0;
    // vRecovery[0] = _recoveredFields[0]->evaluateRecovery(PPR::RECOVERY, 0, x);

    metricData.push_back(vMetric);
    // recoveryData.push_back(vRecovery);
  }

  if(_backmeshOrder > 1) {
    gmsh::view::addModelData(_metricViewTag, 0, bgmToWrite, "NodeData", nodeTagsP1, metricData);
    // gmsh::view::addModelData(recoveryViewTag, 0, bgmToWrite, "NodeData", nodeTags, recoveryData);
  } else {
    gmsh::view::addModelData(_metricViewTag, 0, bgmToWrite, "NodeData", nodeTags, metricData);
    // gmsh::view::addModelData(recoveryViewTag, 0, bgmToWrite, "NodeData", nodeTags, recoveryData);
  }

  // MMG only takes .msh files of version 2.2
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  // Write mesh + metric and mesh + recovery
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);
  // gmsh::view::write(recoveryViewTag, _options.recoveryName);

  // Reset the initial model as active
  gmsh::model::setCurrent(initialMesh);
#else
  UNUSED(nodeTags, coord);
#endif
}

feStatus feMetric::setAnalyticMetric(const std::vector<std::size_t> &nodeTags, const std::vector<double> &coord)
{
  if(_options.analyticMetric == nullptr) {
    return feErrorMsg(FE_STATUS_ERROR, "Cannot set analytic metric because function pointer is null!");
  }

  bool OK = true;
  size_t numVertices = nodeTags.size();

  #if defined(HAVE_OMP)
  #pragma omp parallel shared(OK)
  #endif
  {
    MetricTensor res(1.0);
    feFunctionArguments args;

    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < numVertices; i++)
    {
      if(!OK) continue;

      args.pos[0] = coord[3 * i + 0];
      args.pos[1] = coord[3 * i + 1];
      args.pos[2] = 0.;

      _options.analyticMetric->eval(args, res);

      if(res.determinant() <= 0.) {
        OK = false;
      }

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]].assignMatrixFrom(res);
      }
    }
  }

  if(!OK) {
    return feErrorMsg(FE_STATUS_ERROR, "Determinant is negative or zero for at least one analytic metric tensor!");
  }

  return FE_STATUS_OK;
}

// Compute unscaled optimal metric field minimizing interpolation error
// in Lp norm for P1 elements according to Alauzet & Loseille
feStatus feMetric::computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                                    bool isotropic)
{
  if(_options.useAnalyticDerivatives && _options.secondDerivatives == nullptr) {
    return feErrorMsg(FE_STATUS_ERROR, "Cannot compute analytic derivatives for P1 adaptation because function pointer is null!");
  }

  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double x[2], fxx, fxy, fyx, fyy;
    MetricTensor H, Q1, Q2, Q;
    std::vector<double> D2U_EXACT(4, 0.);
    feFunctionArguments args;

    // Compute bounded absolute value of hessian at vertices (at nodetags)
    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < nodeTags.size(); i++) {
      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

     if(_options.useAnalyticDerivatives) {

        // Evaluate the analytic second derivatives instead
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;
        _options.secondDerivatives->eval(args, D2U_EXACT);

        fxx = D2U_EXACT[0];
        fxy = D2U_EXACT[1];
        fyx = D2U_EXACT[2];
        fyy = D2U_EXACT[3];

      } else {

        // 2nd order derivatives from recoveries
        fxx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 2, x);
        fxy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 3, x);
        fyx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 4, x);
        fyy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 5, x);

      }

      if(_options.targetNorm == Norm::Lp) {

        H(0, 0) = fxx;
        H(0, 1) = (fxy + fyx) / 2.;
        H(1, 0) = (fxy + fyx) / 2.;
        H(1, 1) = fyy;

        if(isotropic)
          Q = H.boundEigenvaluesOfAbsIsotropic(_lambdaMin, _lambdaMax);
        else {
          // Q = H.absoluteValueEigen();
          Q = H.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
        }
      } else {
        // Target norm to minimize is H1. Use Mirebeau's analytical solution.
        double a = 0.5 * fxx;
        double b = 0.5 * (fxy + fyx)/2.;
        double c = 0.5 * fyy;
        // if(fabs(a*c - b*b) < 1e-10) {
        //   feInfo("a*c - b*b = %+-1.4e", a*c - b*b);
        //   exit(-1);
        // }
        double alpha = 100.;
        Q1(0, 0) = 4. * (a*a + b*b);
        Q1(0, 1) = 4. * (a*b + b*c);
        Q1(1, 0) = 4. * (a*b + b*c);
        Q1(1, 1) = 4. * (b*b + c*c);
        // Limit the anisotropy (Mirebeau's thesis, Section 3.4.3 equation 3.74
        Q2 = Q1.limitAnisotropy(alpha);
        Q = Q2.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
      }

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]](0,0) = Q(0,0);
        _metricTensorAtNodetags[nodeTags[i]](0,1) = Q(0,1);
        _metricTensorAtNodetags[nodeTags[i]](1,0) = Q(1,0);
        _metricTensorAtNodetags[nodeTags[i]](1,1) = Q(1,1);
      }
    }
  }

  return FE_STATUS_OK;
}

// Compute unscaled optimal metric field minimizing the error on a goal functional.
// The goal functional is represented through the adjoint problem solution p.
// The optimal metric is the interpolation error metric weighted by the (gradient of the) adjoint
// solution.
feStatus feMetric::computeMetricsGoalOrientedP1(std::vector<std::size_t> &nodeTags,
                                                std::vector<double> &coord)
{
  if(_recoveredFields.size() < 2) {
    return feErrorMsg(
      FE_STATUS_ERROR,
      "Only one recovery was provided to compute the metric field. "
      "The recovered adjoint and its derivatives are also required to compute a metric field for "
      "goal oriented adaptation.");
  }

  double x[2], fxx, fxy, fyx, fyy, p; //, px, py;
  MetricTensor H;

  // Test pour la perturbation singuliere
  // double epsilon = 1e-2;
  // double velocity[2] = {1. - 2. * epsilon, 1. - 2. * epsilon};
  // double a = 2. * (1 - epsilon);

  // Compute bounded absolute value of hessian at vertices (at nodetags)
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    // Evaluate 2nd order derivatives
    fxx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 2, x);
    fxy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 3, x);
    fyx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 4, x);
    fyy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 5, x);

    // Evaluate adjoint and its gradient
    p = _recoveredFields[1]->evaluateRecovery(PPR::RECOVERY, 0, x);
    // px = _recoveredFields[1]->evaluateRecovery(PPR::DERIVATIVE, 0, x);
    // py = _recoveredFields[1]->evaluateRecovery(PPR::DERIVATIVE, 1, x);

    // double dotProd = px * velocity[0] + py * velocity[1];

    H(0, 0) = fxx;
    H(0, 1) = (fxy + fyx) / 2.;
    H(1, 0) = (fxy + fyx) / 2.;
    H(1, 1) = fyy;
    // H *= fabs(dotProd);
    H *= fabs(p);
    _metricTensorAtNodetags[nodeTags[i]] = H.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
  }

  return FE_STATUS_OK;
}

// Compute the anisotropic measure Q for quadratic interpolation error.
// It is the upper bound on the third-order derivatives.
// Analytic computation using Mirebeau's solution.
// Computation in the physical space directly.
feStatus feMetric::computeMetricsP2(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
  if(_options.useAnalyticDerivatives && _options.thirdDerivatives == nullptr) {
    return feErrorMsg(FE_STATUS_ERROR, "Cannot compute analytic derivatives for P1 adaptation because function pointer is null!");
  }

  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_P2) for %d vertices...", nodeTags.size());

  size_t numVertices = nodeTags.size();
  std::map<int, std::vector<double> > &errorCoeffAtVertices = _recoveredFields[0]->getErrorCoefficients();

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double x[2];
    MetricTensor Q;

    // Coefficients of the cubic error polynomial
    size_t numCubicCoeff = 4;
    std::vector<double> errorCoeff(numCubicCoeff, 0.);

    std::vector<double> D3U_EXACT(8, 0.);
    feFunctionArguments args;

// Compute bounded absolute value of upper bound Q at vertices (at nodetags)
#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t i = 0; i < numVertices; i++) {

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      if(_options.useAnalyticDerivatives) {
        // Evaluate the analytic high-order derivatives instead
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;
        _options.thirdDerivatives->eval(args, D3U_EXACT);
        double uxxx = D3U_EXACT[0];
        double uxxy = D3U_EXACT[1];
        double uxyy = D3U_EXACT[3];
        double uyyy = D3U_EXACT[7];
        errorCoeff[0] = uxxx;
        errorCoeff[1] = 3. * uxxy;
        errorCoeff[2] = 3. * uxyy;
        errorCoeff[3] = uyyy;
      } else {

        // Get the coefficients of the homogeneous error polynomial at vertex
        std::vector<double> &errorCoeffAtVertex = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];

        for(size_t j = 0; j < numCubicCoeff; ++j) {
          errorCoeff[j] = errorCoeffAtVertex[j];
        }
      }

      bool success = false;

      if(_options.targetNorm == Norm::Lp) {

        // IMPORTANT: Unlike the other functions, 
        // here the input coefficients are a = uxxx, b = 3*uxxy, c = 3*uxyym d = uyyy
        success = computeAnalyticMetricP2ForLpNorm(errorCoeff, Q, _options.maxAnalyticEllipseDiameter);

        ////////////////////////////////////////////////////////////////
        // Debug plot (must be done here to plot around (x,y)
        // {
        //   double a = errorCoeff[0];
        //   double b = errorCoeff[1];
        //   double c = errorCoeff[2];
        //   double d = errorCoeff[3];
        //   double sc = fmax(fmax(fabs(a), fabs(b)), fmax(fabs(c), fabs(d)));
        //   a /= sc;
        //   b /= sc;
        //   c /= sc;
        //   d /= sc;
        //   std::string name = "debugPlotAnalyticMetric_" + std::to_string(i) + ".pos";
        //   FILE *myfile = fopen(name.data(), "w");
        //   fprintf(myfile, "View \"debugMetric\"{\n");
        //   int N = 200;
        //   double fac = 100.;
        //   double xx, yy, p, v;
        //   for(int iT = 0; iT < N; ++iT) {
        //     xx = cos(2.*M_PI/(double)N * iT); 
        //     yy = sin(2.*M_PI/(double)N * iT);
        //     p = fabs(a*xx*xx*xx + b*xx*xx*yy + c*xx*yy*yy + d*yy*yy*yy);
        //     v = pow(p, -1./3.);
        //     if(fabs(x[0] + v*xx)/fac < 10. && fabs(x[1] + v*yy)/fac < 10.)
        //     fprintf(myfile, "SP(%g,%g,%g){%g};\n", x[0] + (v*xx)/fac, x[1] + (v*yy)/fac, 0., 1.);
        //   }
        //   fprintf(myfile, "};\n"); fclose(myfile);
        // }
        ////////////////////////////////////////////////////////////////

      } else {
        // Upper bound on grad(error)
        success = computeAnalyticMetricP2ForH1semiNorm(errorCoeff, Q);
      }

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        if(success) {
          _metricTensorAtNodetags[nodeTags[i]](0,0) = Q(0,0);
          _metricTensorAtNodetags[nodeTags[i]](0,1) = Q(0,1);
          _metricTensorAtNodetags[nodeTags[i]](1,0) = Q(1,0);
          _metricTensorAtNodetags[nodeTags[i]](1,1) = Q(1,1);
        } else {
          Q(0,0) = _lambdaMin;
          Q(0,1) = 0.;
          Q(1,0) = 0.;
          Q(1,1) = _lambdaMin;
          _metricTensorAtNodetags[nodeTags[i]] = Q.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
          feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e) (vertex %d/%d)", x[0], x[1],
                     i, nodeTags.size());
        }
      }
    }
  }

  return FE_STATUS_OK;
}

feStatus feMetric::computeMetricsP1_referenceSpace(std::vector<std::size_t> &nodeTags, std::vector<double> &coord, int setToIdentity)
{
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_P1 in the REFERENCE SPACE) for %d vertices...", nodeTags.size());

  size_t numVertices = nodeTags.size();
  std::map<int, std::vector<double> > &errorCoeffAtVertices = _recoveredFields[0]->getErrorCoefficients();

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double x[2];
    MetricTensor Qprev, Qtri, Mmud, Mud;
    std::vector<double> D2U_EXACT(4, 0.);
    feFunctionArguments args;

    // Initialize all metrics to identity for the first iteration
    if(setToIdentity) {
      for(size_t i = 0; i < numVertices; i++)
      {
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        {
          _metricTensorAtNodetags[nodeTags[i]](0,0) = 1.;
          _metricTensorAtNodetags[nodeTags[i]](0,1) = 0.;
          _metricTensorAtNodetags[nodeTags[i]](1,0) = 0.;
          _metricTensorAtNodetags[nodeTags[i]](1,1) = 1.;
        }
      }
    }

    // Compute Q_triangle from C_ijk * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko
    // Compute bounded absolute value of upper bound Q at vertices (at nodetags)
    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < numVertices; i++) {

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        Qprev = _metricTensorAtNodetags[nodeTags[i]].copy();
        Mmud = _metricTensorAtNodetags[nodeTags[i]].pow(-0.5);
        Mud  = _metricTensorAtNodetags[nodeTags[i]].sqrt();
      }

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      double uxx, uxy, uyy;

      // Get the coefficients of the homogeneous error polynomial at vertex
      std::vector<double> &errorCoeff = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];
      uxx = errorCoeff[0];
      uxy = errorCoeff[1] / 2.;
      uyy = errorCoeff[2];

      // IMPORTANT: This should be correct and t converges with rate 3 for tanh and a = 10.
      // The error coefficient are multiplied by the corresponding binomial coefficient
      // WHEN EVALUATING THE ERROR POLYNOMIAL, see feMetricTools.cpp.
      // This allows to factorize the homogeneous polynomial.
      // Both methods in "evaluateHomogeneousErrorPolynomial" multiply the error coefficients,
      // so the values in "errorCoeff" should be d^(k+1) u / dx^(k1) dy^(k2), 
      // without any coefficient in front. So for the mixed derivatives, the result
      // stored in errorCoeff must be divided since they were summed in feNewRecovery.

      // Divide the mixed derivatives coefficient by the corresponding binomial coefficient.
      // Coefficients in errorCoeff are the sum of all relevant terms in 
      // the homogeneous polynomial, counting multiple times repeated exponents,
      // e.g. x^2*y appears 3 times as xxy, xyx and yxx.
      // errorCoeff[1] /= 3.;
      // errorCoeff[2] /= 3.;

      if(_options.useAnalyticDerivatives) {
        // Evaluate the analytic high-order derivatives instead
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;
        _options.secondDerivatives->eval(args, D2U_EXACT);
        uxx = D2U_EXACT[0];
        uxy = D2U_EXACT[1];
        uyy = D2U_EXACT[3];
        // errorCoeff[0] = uxx;
        // errorCoeff[1] = 2. * uxy;
        // errorCoeff[2] = uyy;
      }

      // Compute M^(-1/2)
      // MetricTensor Mmud = _metricTensorAtNodetags[nodeTags[i]].pow(-0.5);
      // feInfo("Mmud = %f - %f - %f - %f", Mmud(0,0), Mmud(0,1), Mmud(1,0), Mmud(1,1));

      // Compute C_ijk * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko
      double H[2][2] = {{uxx, uxy}, {uxy, uyy}};
      MetricTensor Hmat(1.);
      Hmat(0,0) = H[0][0];
      Hmat(0,1) = H[0][1];
      Hmat(1,0) = H[1][0];
      Hmat(1,1) = H[1][1];
      MetricTensor absH = Hmat.absoluteValueEigen();

      double res[2][2] = {{0., 0.}, {0., 0.}};
      for(int I = 0; I < 2; ++I) {
        for(int J = 0; J < 2; ++J) {
          for(int m = 0; m < 2; ++m) {
            for(int n = 0; n < 2; ++n) {
              res[I][J] += H[m][n] * Mmud(I,m) * Mmud(J,n);
              // res[I][J] += absH(m,n) * Mmud(I,m) * Mmud(J,n);
            }
          }
        }
      }

      // errorCoeff[0] = res[0][0];
      // errorCoeff[1] = res[0][1] + res[1][0];
      // errorCoeff[3] = res[1][1];
      // feInfo("");
      // feInfo("H     = [%f  %f ; %f  %f]", H[0][0], H[0][1], H[1][0], H[1][1]);
      // feInfo("absH  = [%f  %f ; %f  %f]", absH(0,0), absH(0,1), absH(1,0), absH(1,1));
      // feInfo("Mmud  = [%f  %f ; %f  %f]", Mmud(0,0), Mmud(0,1), Mmud(1,0), Mmud(1,1));
      // feInfo("res   = [%f  %f ; %f  %f]", res[0][0], res[0][1], res[1][0], res[1][1]);

      // computeAnalyticMetricP2ForLpNorm(errorCoeff, Qtri, _options.maxAnalyticEllipseDiameter);
      MetricTensor mat(1.);
      mat(0,0) = res[0][0];
      mat(0,1) = res[0][1];
      mat(1,0) = res[1][0];
      mat(1,1) = res[1][1];
      // feInfo("mat   = [%f  %f ; %f  %f]", mat(0,0), mat(0,1), mat(1,0), mat(1,1));

      Qtri = mat.absoluteValueEigen();
      // feInfo("Qtri  = [%f  %f ; %f  %f]", Qtri(0,0), Qtri(0,1), Qtri(1,0), Qtri(1,1));

      MetricTensor Q = Mud * Qtri * Mud;

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]](0,0) = Q(0,0);
        _metricTensorAtNodetags[nodeTags[i]](0,1) = Q(0,1);
        _metricTensorAtNodetags[nodeTags[i]](1,0) = Q(1,0);
        _metricTensorAtNodetags[nodeTags[i]](1,1) = Q(1,1);
      }
    }
  }
  return FE_STATUS_OK;
}

// Compute the anisotropic measure Q for quadratic interpolation error.
// Same as computeMetricsP2, but done in the reference space iteratively.
feStatus feMetric::computeMetricsP2_referenceSpace(std::vector<std::size_t> &nodeTags, std::vector<double> &coord, int setToIdentity)
{
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_P2 in the REFERENCE SPACE) for %d vertices...", nodeTags.size());

  size_t numVertices = nodeTags.size();
  std::map<int, std::vector<double> > &errorCoeffAtVertices = _recoveredFields[0]->getErrorCoefficients();

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double x[2];
    MetricTensor Qprev, Qtri, Mmud, Mud, foo;
    std::vector<double> D3U_EXACT(8, 0.);
    feFunctionArguments args;

    // Initialize all metrics to identity for the first iteration
    if(setToIdentity) {
      for(size_t i = 0; i < numVertices; i++)
      {
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        {
          _metricTensorAtNodetags[nodeTags[i]](0,0) = 1.;
          _metricTensorAtNodetags[nodeTags[i]](0,1) = 0.;
          _metricTensorAtNodetags[nodeTags[i]](1,0) = 0.;
          _metricTensorAtNodetags[nodeTags[i]](1,1) = 1.;
        }
      }
    }

    // Compute Q_triangle from C_ijk * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko
    // Compute bounded absolute value of upper bound Q at vertices (at nodetags)
    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < numVertices; i++) {

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        Qprev = _metricTensorAtNodetags[nodeTags[i]].copy();
        Mmud = _metricTensorAtNodetags[nodeTags[i]].pow(-0.5);
        Mud  = _metricTensorAtNodetags[nodeTags[i]].sqrt();
        foo  = _metricTensorAtNodetags[nodeTags[i]].pow(0.75);
      }

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      double uxxx, uxxy, uxyy, uyyy;

      // Get the coefficients of the homogeneous error polynomial at vertex
      std::vector<double> &errorCoeff = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];
      uxxx = errorCoeff[0];
      uxxy = errorCoeff[1] / 3.;
      uxyy = errorCoeff[2] / 3.;
      uyyy = errorCoeff[3];

      // IMPORTANT: This should be correct and t converges with rate 3 for tanh and a = 10.
      // The error coefficient are multiplied by the corresponding binomial coefficient
      // WHEN EVALUATING THE ERROR POLYNOMIAL, see feMetricTools.cpp.
      // This allows to factorize the homogeneous polynomial.
      // Both methods in "evaluateHomogeneousErrorPolynomial" multiply the error coefficients,
      // so the values in "errorCoeff" should be d^(k+1) u / dx^(k1) dy^(k2), 
      // without any coefficient in front. So for the mixed derivatives, the result
      // stored in errorCoeff must be divided since they were summed in feNewRecovery.

      // Divide the mixed derivatives coefficient by the corresponding binomial coefficient.
      // Coefficients in errorCoeff are the sum of all relevant terms in 
      // the homogeneous polynomial, counting multiple times repeated exponents,
      // e.g. x^2*y appears 3 times as xxy, xyx and yxx.
      // errorCoeff[1] /= 3.;
      // errorCoeff[2] /= 3.;

      if(_options.useAnalyticDerivatives) {
        // Evaluate the analytic high-order derivatives instead
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;
        _options.thirdDerivatives->eval(args, D3U_EXACT);
        uxxx = D3U_EXACT[0];
        uxxy = D3U_EXACT[1];
        uxyy = D3U_EXACT[3];
        uyyy = D3U_EXACT[7];
        errorCoeff[0] = uxxx;
        errorCoeff[1] = 3. * uxxy;
        errorCoeff[2] = 3. * uxyy;
        errorCoeff[3] = uyyy;
      }

      // // feInfo("Coeff 1 = %f - %f - %f - %f", errorCoeff[0], errorCoeff[1], errorCoeff[2], errorCoeff[3]);

      // Compute M^(-1/2)
      // MetricTensor Mmud = _metricTensorAtNodetags[nodeTags[i]].pow(-0.5);
      // feInfo("Mmud = %f - %f - %f - %f", Mmud(0,0), Mmud(0,1), Mmud(1,0), Mmud(1,1));

      // Compute C_ijk * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko
      double C[2][2][2] = {{{uxxx, uxxy}, {uxxy, uxyy}}, {{uxxy, uxyy}, {uxyy, uyyy}}};
      double res[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
      for(int I = 0; I < 2; ++I) {
        for(int J = 0; J < 2; ++J) {
          for(int K = 0; K < 2; ++K) {
            for(int m = 0; m < 2; ++m) {
              for(int n = 0; n < 2; ++n) {
                for(int o = 0; o < 2; ++o) {
                  res[I][J][K] += C[m][n][o] * Mmud(I,m) * Mmud(J,n) * Mmud(K,o);
                }
              }
            }
          }
        }
      }

      errorCoeff[0] = res[0][0][0];
      errorCoeff[1] = res[0][0][1] + res[0][1][0] + res[1][0][0];
      errorCoeff[2] = res[0][1][1] + res[1][0][1] + res[1][1][0];
      errorCoeff[3] = res[1][1][1];

      // feInfo("Coeff 2 = %f - %f - %f - %f", errorCoeff[0], errorCoeff[1], errorCoeff[2], errorCoeff[3]);

      computeAnalyticMetricP2ForLpNorm(errorCoeff, Qtri, _options.maxAnalyticEllipseDiameter);

      // MetricTensor tmp = Qtri * Mud;
      // MetricTensor Q = Mud * tmp;
      MetricTensor Q = Mud * Qtri * Mud;

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]](0,0) = Q(0,0);
        _metricTensorAtNodetags[nodeTags[i]](0,1) = Q(0,1);
        _metricTensorAtNodetags[nodeTags[i]](1,0) = Q(1,0);
        _metricTensorAtNodetags[nodeTags[i]](1,1) = Q(1,1);
      }
    }
  }
  return FE_STATUS_OK;
}

#if defined(HAVE_SOPLEX)
void setUpLinearProblem(linearProblem &myLP, feMetricOptions &options, int nTheta,
                        bool reset = false)
{
  if(reset) {
    myLP.problem.clearLPReal();
    myLP.lprowset.clear();
  }
  // Set up the optimization problem structure
  int numConstraints = 4 * nTheta;
  size_t size = 2 * numConstraints; // 2 coordonnees * 4 quadrants * nTheta angles/quadrant
  myLP.lvl1.resize(size, 0.);
  myLP.constraints.resize(size, 0.);
  myLP.numConstraints = numConstraints;
  myLP.uniformErrorCurve = options.logSimplexOptions.uniformErrorCurve;
  myLP.numLoopsUniformErrorCurve = options.logSimplexOptions.numLoopsUniformErrorCurve;

  // Set generic LP solver data
  myLP.problem.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
  myLP.problem.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MINIMIZE);
  DSVector dummycol(0);
  double infty = 1e+100;
  myLP.problem.addColReal(LPCol(1.0, dummycol, infty, -infty));
  myLP.problem.addColReal(LPCol(0.0, dummycol, infty, -infty));
  myLP.problem.addColReal(LPCol(1.0, dummycol, infty, -infty));

  myLP.row.setMax(3);
  myLP.prim.reDim(3);
  // myLP.lprow = LPRow(3);
  // myLP.lprow.setRhs(infty);
  myLP.lprowset = LPRowSet(numConstraints, 3);

  for(int i = 0; i < numConstraints; ++i) {
    DSVector row(3);
    row.add(0, 1.0);
    row.add(1, 1.0);
    row.add(2, 1.0);
    myLP.problem.addRowReal(LPRow(1.0, row, infty));
  }

  // Used to compute the roots of the degree 6 error polynomial for curved adaptation
  myLP.lvl1PolynomialCoeffs_minus = Eigen::VectorXd::Zero(7);
  myLP.lvl1PolynomialCoeffs_plus = Eigen::VectorXd::Zero(7);
}
#endif

#if defined(TEST_ANALYTIC_METRIC)
static MetricTensor metricAparicioEstrems(double x, double y)
{
  double A = 10.;
  double c = sqrt(A * A + 4. * M_PI * M_PI);
  double HMIN = 0.1;
  double b = 10.;

  // Size from the paper
  // double h = HMIN + b * fabs(y);
  // double dhdx = 0.;
  // double dhdy = b * sgn(y);

  double h = HMIN + b * fabs((A * y - cos(2. * M_PI * x)) / c);
  // double dhdx = b * (2. * pi * sin(2. * pi * x) * sgn(A * y - cos(2. * pi * x))) / c;
  // double dhdy = b * (A * sgn(A * y - cos(2. * pi * x))) / c;

  MetricTensor D(0.);
  D(0, 0) = 1.;
  D(1, 1) = 1. / (h * h);

  // Gradient of (x, (10y - cos(2*pi*x)) / c)
  // FIXME: gradPhi should be a SquareMatrix, not a MetricTensor
  MetricTensor gradPhi(0.);
  gradPhi(0, 0) = 1.;
  gradPhi(0, 1) = 0.;
  gradPhi(1, 0) = 2. * M_PI * sin(2. * M_PI * x) / c;
  gradPhi(1, 1) = A / c;

  return gradPhi.transpose() * D * gradPhi;
}
#endif

feStatus feMetric::computeMetricsP2_forGraphSurface(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
    if(_options.useAnalyticDerivatives && _options.thirdDerivatives == nullptr) {
    return feErrorMsg(FE_STATUS_ERROR, "Cannot compute analytic derivatives for adaptation because function pointer is null!");
  }

  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_P2_GRAPH_SURFACE) for %d vertices...", nodeTags.size());

  size_t numVertices = nodeTags.size();
  // std::map<int, std::vector<double> > &errorCoeffAtVertices = _recoveredFields[0]->getErrorCoefficients();

  FILE *myFile = fopen("graphSurface.pos", "w");
  fprintf(myFile, "View \" graphSurface \"{ \n");

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double x[2];
    MetricTensor Q;
    std::vector<double> DU_EXACT(2, 0.), D2U_EXACT(4, 0.), D3U_EXACT(8, 0.);
    feFunctionArguments args;

// Compute bounded absolute value of upper bound Q at vertices (at nodetags)
#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t i = 0; i < numVertices; i++) {

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      if(_options.useAnalyticDerivatives) {
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;

        _options.firstDerivatives->eval(args, DU_EXACT);
        _options.secondDerivatives->eval(args, D2U_EXACT);
        _options.thirdDerivatives->eval(args, D3U_EXACT);
      } else {
        feInfo("Implement GRAPH_SURFACE metric for recovered derivatives.");
        exit(-1);
      }

      computeGraphSurfaceMetric(x, DU_EXACT.data(), D2U_EXACT.data(), D3U_EXACT.data(), Q, myFile);

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]](0,0) = Q(0,0);
        _metricTensorAtNodetags[nodeTags[i]](0,1) = Q(0,1);
        _metricTensorAtNodetags[nodeTags[i]](1,0) = Q(1,0);
        _metricTensorAtNodetags[nodeTags[i]](1,1) = Q(1,1);
      }
    }
  }

  fprintf(myFile, "};"); fclose(myFile);

  return FE_STATUS_OK;
}

// Compute scaled optimal metric field for Pn elements
// according to Coulaud & Loseille using the log-simplex method
feStatus feMetric::computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                                    bool isotropic)
{
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_PN) for %d vertices...", nodeTags.size());

  if(_options.polynomialDegree > 4) {
    return feErrorMsg(FE_STATUS_ERROR,
                      "Modify mixed derivative error coefficients for polynomials of order %d",
                      _options.polynomialDegree);
  }

#if defined(ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE)
  // Just assign analytic metric and exit
  // To test convergence of metric interpolation
  for(size_t i = 0; i < nodeTags.size(); i++) {
    double X = coord[3 * i + 0];
    double Y = coord[3 * i + 1];
    _metricTensorAtNodetags[nodeTags[i]] = analyticMetric(X, Y);
  }
  return FE_STATUS_OK;
#endif

#if defined(TEST_ANALYTIC_METRIC)
  // Just assign analytic metric and exit
  // To test distortion minimization
  for(size_t i = 0; i < nodeTags.size(); i++) {
    double X = coord[3 * i + 0];
    double Y = coord[3 * i + 1];
    _metricTensorAtNodetags[nodeTags[i]] = metricAparicioEstrems(X, Y);
  }
  return FE_STATUS_OK;
#endif

#if defined(HAVE_SOPLEX)

  std::map<int, std::vector<double> > &errorCoeffAtVertices =
    _recoveredFields[0]->getErrorCoefficients();

  bool OK = true;

  size_t numVertices = nodeTags.size();

  double n = 2.; // Space dimension, fixed for now (-:

#if defined(HAVE_OMP)
#pragma omp parallel shared(OK)
#endif
  {
    double x[2];
    MetricTensor Q;
    
    // Coefficients of the error polynomial
    std::vector<double> errorCoeff(_options.polynomialDegree + 2, 0.);

    // Coefficients of the squared norm of the gradient of the error polynomial
    // Gradient of error polynomial SQUARED is homogeneous of degree 2*k, and has 2*k+1 coefficients
    std::vector<double> coeffGradErrorPolynomial(n * _options.polynomialDegree + 1, 0.);

    std::vector<double> D2U_EXACT(4, 0.), D3U_EXACT(8, 0.), D4U_EXACT(16, 0.), D5U_EXACT(32, 0.);
    feFunctionArguments args;

    int numIter;
    int nTheta = _options.logSimplexOptions.nThetaPerQuadrant;
    int maxIter = _options.logSimplexOptions.maxIter;
    double tol = _options.logSimplexOptions.tol;

    linearProblem myLP;
    setUpLinearProblem(myLP, _options, nTheta);

    // Compute bounded absolute value of upper bound Q at vertices (at nodetags)
    #if defined(HAVE_OMP)
    #pragma omp for schedule(dynamic, 1)
    #endif
    for(size_t i = 0; i < numVertices; i++) {
      if(!OK) continue;

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      if(_options.useAnalyticDerivatives) {
        args.pos[0] = x[0];
        args.pos[1] = x[1];
        args.pos[2] = 0.;
        // Evaluate the analytic high-order derivatives instead
        if(_options.polynomialDegree == 1) {
          
          _options.secondDerivatives->eval(args, D2U_EXACT);

          double uxx = D2U_EXACT[0];
          double uxy = D2U_EXACT[1];
          double uyy = D2U_EXACT[3];

          errorCoeff[0] = uxx;
          errorCoeff[1] = uxy;
          errorCoeff[2] = uyy;

        } else if(_options.polynomialDegree == 2) {
          
          _options.thirdDerivatives->eval(args, D3U_EXACT);

          double uxxx = D3U_EXACT[0];
          double uxxy = D3U_EXACT[1];
          double uxyy = D3U_EXACT[3];
          double uyyy = D3U_EXACT[7];

          errorCoeff[0] = uxxx;
          errorCoeff[1] = uxxy;
          errorCoeff[2] = uxyy;
          errorCoeff[3] = uyyy;

        } else if(_options.polynomialDegree == 3) {

          _options.fourthDerivatives->eval(args, D4U_EXACT);

          errorCoeff[0] = D4U_EXACT[0];
          errorCoeff[1] = D4U_EXACT[1];
          errorCoeff[2] = D4U_EXACT[3];
          errorCoeff[3] = D4U_EXACT[7];
          errorCoeff[4] = D4U_EXACT[15];

        } else if(_options.polynomialDegree == 4) {

          _options.fifthDerivatives->eval(args, D5U_EXACT);

          errorCoeff[0] = D5U_EXACT[0];
          errorCoeff[1] = D5U_EXACT[1];
          errorCoeff[2] = D5U_EXACT[3];
          errorCoeff[3] = D5U_EXACT[7];
          errorCoeff[4] = D5U_EXACT[15];
          errorCoeff[5] = D5U_EXACT[31];

        } else {
          feInfo("TEMPORARY ERROR MESSAGE: CANNOT COMPUTE ANALYTIC DERIVATIVES FOR THIS POLYNOMIAL DEGREE");
          exit(-1);
        }
      } else {

        // Get the coefficients of the homogeneous error polynomial at vertex
        std::vector<double> &errorCoeffAtVertex = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];

        for(size_t j = 0; j < _options.polynomialDegree + 1; ++j) {
          errorCoeff[j] = errorCoeffAtVertex[j];
        }

        // IMPORTANT: This should be correct and t converges with rate 3 for tanh and a = 10.
        // The error coefficient are multiplied by the corresponding binomial coefficient
        // WHEN EVALUATING THE ERROR POLYNOMIAL, see feMetricTools.cpp.
        // This allows to factorize the homogeneous polynomial.
        // Both methods in "evaluateHomogeneousErrorPolynomial" multiply the error coefficients,
        // so the values in "errorCoeff" should be d^(k+1) u / dx^(k1) dy^(k2), 
        // without any coefficient in front. So for the mixed derivatives, the result
        // stored in errorCoeff must be divided since they were summed in feNewRecovery.

        // Divide the mixed derivatives coefficient by the corresponding binomial coefficient.
        // Coefficients in errorCoeff are the sum of all relevant terms in 
        // the homogeneous polynomial, counting multiple times repeated exponents,
        // e.g. x^2*y appears 3 times as xxy, xyx and yxx.
        if(_options.polynomialDegree == 1) {
          errorCoeff[1] /= 2.;
        } else if(_options.polynomialDegree == 2) {
          errorCoeff[1] /= 3.;
          errorCoeff[2] /= 3.;
        } else if(_options.polynomialDegree == 3) {
          errorCoeff[1] /= 4.;
          errorCoeff[2] /= 6.;
          errorCoeff[3] /= 4.;
        } else if(_options.polynomialDegree == 4) {
          errorCoeff[1] /= 5.;
          errorCoeff[2] /= 10.;
          errorCoeff[3] /= 10.;
          errorCoeff[4] /= 5.;
        } else {
          feInfo("TEMPORARY ERROR MESSAGE: IMPLEMENT ERROR COEFFICIENTS FOR THIS POLYNOMIAL DEGREE");
          exit(-1);
        }
      }

      bool res = false;

      // This is the degree-1 of the ERROR POLYNOMIAL
      // I do this degree + 1 when evaluating the polynomial in the logsimplex 
      int polynomialDegreeMinusOne = _options.polynomialDegree;
      bool applyBinomialCoefficients = true;

      if(_options.targetNorm == Norm::Lp) {
        // Optimal upper bound for |pi(z)| <= 1
        res = computeMetricLogSimplexStraight(x, errorCoeff, polynomialDegreeMinusOne, nTheta,
                                                 maxIter, tol, applyBinomialCoefficients, Q, numIter, myLP);
      } else {
        // Optimal upper bound for |grad(pi)(z)| <= 1

        // The coefficients passed to the logsimplex function
        // are the exact coefficients of the gradient of the 
        // error polynomial. No need to multiply by binomial coefficients
        // when evaluating the polynomial.
        applyBinomialCoefficients = false;

        // |grad pi|^2 = dx_pi^2 + dy_pi^2 is of degree dim*k
        polynomialDegreeMinusOne = n * _options.polynomialDegree - 1; 

        if(_options.polynomialDegree == 1) {

          // The coefficients of the SQUARED gradient of 
          // pi(x,y) = ax^2 + bxy + cy^2
          // The coefficients above include the multiple, e.g.:
          // a =     d_xx f
          // b = 2 * d_xy f
          // c =     d_yy f

          double a =      errorCoeff[0];
          double b = 2. * errorCoeff[1];
          double c =      errorCoeff[2];

          coeffGradErrorPolynomial[0] = 4 * a*a +     b*b; // x^2
          coeffGradErrorPolynomial[1] = 4 * a*b + 4 * b*c; // x y
          coeffGradErrorPolynomial[2] =     b*b + 4 * c*c; // y^2

        } else if(_options.polynomialDegree == 3) {

          // The coefficients of the SQUARED gradient of 
          // pi(x,y) = ax^4 + bx^3y + cx^2y^2 + dxy^3 + ey^4
          // The coefficients above include the multiple, e.g.:
          // a =     d_xxxx f
          // b = 4 * d_xxxy f
          // c = 6 * d_xxyy f
          // d = 4 * d_xyyy f
          // e =     d_yyyy f

          double a =      errorCoeff[0];
          double b = 4. * errorCoeff[1];
          double c = 6. * errorCoeff[2];
          double d = 4. * errorCoeff[3];
          double e =      errorCoeff[4];

          coeffGradErrorPolynomial[0] = 16 * a*a +      b*b; // x^6
          coeffGradErrorPolynomial[1] = 24 * a*b +  4 * b*c; // x^5 y
          coeffGradErrorPolynomial[2] =  9 * b*b +  6 * b*d + 4 * c*c + 16 * a*c; // x^4 y^2
          coeffGradErrorPolynomial[3] =  8 * a*d + 12 * b*c + 8 * b*e + 12 * c*d; // x^3 y^3
          coeffGradErrorPolynomial[4] =  9 * d*d +  6 * b*d + 4 * c*c + 16 * e*c; // x^2 y^4
          coeffGradErrorPolynomial[5] = 24 * d*e +  4 * c*d; // x y^5
          coeffGradErrorPolynomial[6] = 16 * e*e +      d*d; // y^6

        } else if(_options.polynomialDegree == 4) {

          // The coefficients of the SQUARED gradient of 
          // pi(x,y) = ax^5 + bx^4y + cx^3y^2 + dx^2y^3 + exy^4 + fy^5
          // The coefficients above include the multiple, e.g.:
          // a =      d_xxxxx f
          // b =  5 * d_xxxxy f
          // c = 10 * d_xxxyy f
          // d = 10 * d_xxyyy f
          // e =  5 * d_xyyyy f
          // f =      d_yyyyy f

          double a =       errorCoeff[0];
          double b =  5. * errorCoeff[1];
          double c = 10. * errorCoeff[2];
          double d = 10. * errorCoeff[3];
          double e =  5. * errorCoeff[4];
          double f =       errorCoeff[5];

          coeffGradErrorPolynomial[0] = 25 * a*a +      b*b; // x^8
          coeffGradErrorPolynomial[1] = 40 * a*b +  4 * b*c; // x^7 y
          coeffGradErrorPolynomial[2] = 16 * b*b +  6 * b*d + 30 * a*c + 4 * c*c; // x^6 y^2
          coeffGradErrorPolynomial[3] = 20 * a*d + 24 * b*c + 8 * b*e + 12 * c*d; // x^5 y^3
          coeffGradErrorPolynomial[4] =  9 * (c*c + d*d) + 16 * (e*c + b*d) + 10 * (a*e + b*f); // x^4 y^4
          coeffGradErrorPolynomial[5] = 20 * c*f + 24 * d*e + 8 * b*e + 12 * c*d; // x^3 y^5
          coeffGradErrorPolynomial[6] = 16 * e*e +  6 * c*e + 30 * f*d + 4 * d*d; // x^2 y^6
          coeffGradErrorPolynomial[7] = 40 * e*f +  4 * d*e; // x y^7
          coeffGradErrorPolynomial[8] = 25 * f*f +      e*e; // y^8

        } else {
          OK = false;
          feErrorMsg(FE_STATUS_ERROR, "Implement optimal metric for H1 seminorm for P4+");
        }

        res = computeMetricLogSimplexStraight(x, coeffGradErrorPolynomial, polynomialDegreeMinusOne, nTheta,
                                                 maxIter, tol, applyBinomialCoefficients, Q, numIter, myLP);

        // double aa = errorCoeff[0];
        // double bb = errorCoeff[1];
        // double cc = errorCoeff[2];

        // double alpha = 100.;
        // MetricTensor Q1(1.0);
        // Q1(0, 0) = 4. * (aa*aa + bb*bb);
        // Q1(0, 1) = 4. * (aa*bb + bb*cc);
        // Q1(1, 0) = 4. * (aa*bb + bb*cc);
        // Q1(1, 1) = 4. * (bb*bb + cc*cc);

        // // Limit the anisotropy (Mirebeau's thesis, Section 3.4.3 equation 3.74
        // MetricTensor Q2 = Q1.limitAnisotropy(alpha);
        // // Q = H.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);

        // MetricTensor diff1 = Q - Q1;
        // MetricTensor diff2 = Q - Q2;


        // feInfo("LS solution  : %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", Q(0,0), Q(0,1), Q(1,0), Q(1,1));
        // feInfo("Ex solution 1: %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - det = %+-1.3e", Q1(0,0), Q1(0,1), Q1(1,0), Q1(1,1), 100. * diff1.frobeniusNorm()/Q.frobeniusNorm(), aa*cc - bb*bb);
        // feInfo("Ex solution 2: %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - det = %+-1.3e", Q2(0,0), Q2(0,1), Q2(1,0), Q2(1,1), 100. * diff2.frobeniusNorm()/Q.frobeniusNorm(), aa*cc - bb*bb);

        // FILE *myfile = fopen("debug_H1.pos", "w");
        // fprintf(myfile, "View\"debug_H1\"{\n");
        // drawSingleEllipse(myfile, x, Q, 1, 30);
        // drawSingleEllipse(myfile, x, Q1, 1, 30);

        // aa = errorCoeff[0];
        // bb = 2. * errorCoeff[1];
        // cc = errorCoeff[2];

        // MetricTensor Q3(1.0);
        // Q3(0, 0) = 4. * (aa*aa + bb*bb);
        // Q3(0, 1) = 4. * (aa*bb + bb*cc);
        // Q3(1, 0) = 4. * (aa*bb + bb*cc);
        // Q3(1, 1) = 4. * (bb*bb + cc*cc);

        // drawSingleEllipse(myfile, x, Q3, 1, 30);
        // fprintf(myfile, "};\n"); fclose(myfile);

        // MetricTensor diff3 = Q - Q3;
        // feInfo("Ex solution 3: %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - det = %+-1.3e", Q3(0,0), Q3(0,1), Q3(1,0), Q3(1,1), 100. * diff3.frobeniusNorm()/Q.frobeniusNorm(), aa*cc - bb*bb);

        // MetricTensor Q4 = Q3.sqrt();
        // MetricTensor diff4 = Q - Q4;
        // feInfo("Ex solution 4: %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - det = %+-1.3e", Q4(0,0), Q4(0,1), Q4(1,0), Q4(1,1), 100. * diff4.frobeniusNorm()/Q.frobeniusNorm(), aa*cc - bb*bb);
      }

      // If metric computation failed, try again by doubling the number of constraints points in
      // each quadrant (x8 total)
      int retry = 0, maxTry = 5;
      while(!res && retry < maxTry) {
        nTheta *= 2;
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
                   "Could not compute metric. Trying again with %d constraints.", nTheta);
        setUpLinearProblem(myLP, _options, nTheta, true);
        res = computeMetricLogSimplexStraight(x, errorCoeff, polynomialDegreeMinusOne, nTheta,
                                              maxIter, tol, applyBinomialCoefficients, Q, numIter, myLP);
        if(res) feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Success");
        retry++;
      }

      // Restore initial number of constraints
      nTheta = _options.logSimplexOptions.nThetaPerQuadrant;

      if(res) {
        // Bound the eigenvalues of Q
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        {
          if(isotropic) {
            _metricTensorAtNodetags[nodeTags[i]] =
              Q.boundEigenvaluesOfAbsIsotropic(_lambdaMin, _lambdaMax);
          } else
            _metricTensorAtNodetags[nodeTags[i]] = Q.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
        }

        // feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
        //        "Computed metric in %2d iterations - vertex %6d/%6d",
        //        numIter, ++cnter, numVertices);
      } else {
        // OK = false;
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        Q(0,0) = _lambdaMin;
        Q(0,1) = 0.;
        Q(1,0) = 0.;
        Q(1,1) = _lambdaMin;
        _metricTensorAtNodetags[nodeTags[i]] = Q.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
        feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e) (vertex %d/%d)", x[0], x[1],
                   i, nodeTags.size());
      }
    }
  }

  if(!OK) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not compute at least one metric tensor");
  }

  return FE_STATUS_OK;
#else
  UNUSED(coord, isotropic);
  return feErrorMsg(
    FE_STATUS_ERROR,
    "SoPlex (optimization library) is required to compute metric tensors for Pn adaptation");
#endif
}

 void feMetric::computeSizeField(const double pos[2], const int vertex, const double directionGrad[2],
  double &hGrad, double &hIso, const double *du, const double *d2u, const double *d3u)
{
  UNUSED(pos, vertex);
  switch(_options.polynomialDegree) {
    case 1:
      {
        // double derivativeAlongIsoline = fabs(secondDerivativesAlongCurve(pos, vertex, directionGrad, 0));
        // double derivativeAlongGradient = fabs(secondDerivativesAlongCurve(pos, vertex, directionGrad, 1));
        // hIso  = sqrt(2.0 * _options.eTargetError / fmax(derivativeAlongIsoline, 1e-10));
        // hGrad = sqrt(2.0 * _options.eTargetError / fmax(derivativeAlongGradient, 1e-10));

        // Test for f(x,y) = x^2+y^2 only
        // hGrad = sqrt(_options.eTargetError);
        // // hIso = sqrt(2. * _options.eTargetError * sqrt(pos[0]*pos[0]+pos[1]*pos[1]));
        // hIso = sqrt(sqrt(6.)) * sqrt(_options.eTargetError * sqrt(pos[0]*pos[0]+pos[1]*pos[1]));

        // hIso  = solveSizePolynomialLinear(pos, vertex, directionGrad, 0, pow(_options.eTargetError, 2) );
        // hGrad = solveSizePolynomialLinear(pos, vertex, directionGrad, 1,     _options.eTargetError     );

        // hIso  = solveSizePolynomialLinear(pos, vertex, directionGrad, 0, _options.eTargetError, du, d2u);
        // hGrad = solveSizePolynomialLinear(pos, vertex, directionGrad, 1, _options.eTargetError, du, d2u);

        solveSizePolynomialLinear(_options.eTargetError, directionGrad, du, d2u, hIso, hGrad);

        // if(fabs(hIso_ini - hIso)/fabs(hIso_ini) > 1e-4 || fabs(hGrad_ini - hGrad)/fabs(hGrad_ini) > 1e-4) {
        //   feInfo("At (%+-1.4e,%+-1.4e) : %+-1.4e vs %+-1.4e -- %+-1.4e vs %+-1.4e", pos[0], pos[1], hIso_ini, hIso, hGrad_ini, hGrad);
        //   feErrorMsg(FE_STATUS_ERROR, "Human is dead. Mismatch.");
        //   exit(-1);
        // }

      }
      return;
    case 2:
      {
        // double derivativeAlongGradient = fabs(thirdDerivativesAlongCurve(pos, vertex, directionGrad, 0));
        // double derivativeAlongIsoline = fabs(thirdDerivativesAlongCurve(pos, vertex, directionGrad, 1));
        // hIso = pow(6.0 * _options.eTargetError / derivativeAlongIsoline, 1./3.);
        // hGrad = pow(6.0 * _options.eTargetError / derivativeAlongGradient, 1./3.);

        solveSizePolynomialQuadratic(_options.eTargetError, directionGrad, du, d2u, d3u, hIso, hGrad);

      }
      return;
    default:
      feErrorMsg(FE_STATUS_ERROR, "Cannot compute size field for interpolation of order > 2. Exiting.");
      exit(-1);
  }
}

thread_local std::vector<double> DU(2, 0.), D2U(4, 0.), D3U(8, 0.), D4U(16, 0.), D5U(32, 0.);
thread_local feFunctionArguments ARGS;;

void feMetric::getDerivativesAtVertex(const int derivativesOrder, const double *x, const size_t vertex,
  feNewRecovery *recoveredField, double *du, double *d2u, double *d3u, double *d4u, double *d5u)
{
  if(_options.useAnalyticDerivatives) {

    // Evaluate the analytic gradient
    ARGS.pos[0] = x[0];
    ARGS.pos[1] = x[1];
    ARGS.pos[2] = 0.;

    if(derivativesOrder > 0) {
      _options.firstDerivatives->eval(ARGS, DU);
      du[0] = DU[0];
      du[1] = DU[1];
    }

    if(derivativesOrder > 1) {
      _options.secondDerivatives->eval(ARGS, D2U);
      for(int i = 0; i < 4; ++i)
        d2u[i] = D2U[i];
    }

    if(derivativesOrder > 2) {
      _options.thirdDerivatives->eval(ARGS, D3U);
      for(int i = 0; i < 8; ++i)
        d3u[i] = D3U[i];
    }

    if(derivativesOrder > 3) {
      _options.fourthDerivatives->eval(ARGS, D4U);
      for(int i = 0; i < 16; ++i)
        d4u[i] = D4U[i];
    }

    if(derivativesOrder > 4) {
      _options.fifthDerivatives->eval(ARGS, D5U);
      for(int i = 0; i < 32; ++i)
        d5u[i] = D5U[i];
    }

  } else {

    // Get recovered derivatives
    if(derivativesOrder > 0) {
      du[0] = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
      du[1] = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);
    }

    if(derivativesOrder > 1) {
      d2u[0] = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
      d2u[1] = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) + 
                recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) / 2.;
      d2u[2] = d2u[1];
      d2u[3] = recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);
    }

    if(derivativesOrder > 2) {
      double fxxx =  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
      double fxxy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex) +
                     recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex) +
                     recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex)) / 3.;
      double fxyy = (recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex) +
                     recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex) +
                     recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex)) / 3.;
      double fyyy =  recoveredField->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);
      d3u[0] = fxxx;
      d3u[1] = fxxy;
      d3u[2] = fxxy;
      d3u[3] = fxyy;
      d3u[4] = fxxy;
      d3u[5] = fxyy;
      d3u[6] = fxyy;
      d3u[7] = fyyy;
    }

    if(derivativesOrder > 3) {
      feInfo("Add derivatives of order 4 here");
      exit(-1);
      // res[0] = fxxxx;
      // res[1] = fxxxy;
      // res[2] = fxxxy;
      // res[3] = fxxyy;
      // res[4] = fxxxy;
      // res[5] = fxxyy;
      // res[6] = fxxyy;
      // res[7] = fxyyy;
      // res[8] = fxxxy;
      // res[9] = fxxyy;
      // res[10] = fxxyy;
      // res[11] = fxyyy;
      // res[12] = fxxyy;
      // res[13] = fxyyy;
      // res[14] = fxyyy;
      // res[15] = fyyyy;
    } 
  }
}


feStatus feMetric::computeMetricsCurvedIsolines(const std::vector<std::size_t> &nodeTags,
                                                const std::vector<double> &coord)
{
  size_t numVertices = nodeTags.size();

  FILE *directions = nullptr, *dirSmoothed = nullptr;
  if(_options.debug) {
    directions = fopen("directions.pos", "w");
    dirSmoothed = fopen("directionsSmoothed.pos", "w");
    fprintf(directions, "View \" directions \"{\n");
    fprintf(dirSmoothed, "View \" directionsSmoothed \"{\n");
  }

  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  // Compute the direction field
  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double x[2], directionsGrad[2], C, S, du[2], d2u[4], d3u[8];

    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < nodeTags.size(); i++) {

      // Parallel read in shared map: should be thread-safe?
      int tag = _nodeTag2sequentialTag[nodeTags[i]];

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      getDerivativesAtVertex(1, x, tag, _recoveredFields[0], du);

      computeDirectionFieldFromGradient(x, tag, du, 1e-8, directionsGrad);
      C = directionsGrad[0];
      S = directionsGrad[1];

      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        COS[nodeTags[i]] = C;
        SIN[nodeTags[i]] = S;
      }
      
      if(_options.debug) {
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        {
          fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
          fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -1.5*S, 1.5*C, 0.);
        }
      }
    }

    if(_options.smoothDirectionField) {
      // Smooth the direction field
      int smoothMaxIter = 100;
      double smoothTol = 0.95;
      smoothDirections(COS, SIN,  smoothMaxIter, smoothTol);
    }

    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < numVertices; i++) {

      int tag = _nodeTag2sequentialTag[nodeTags[i]];

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      C = COS[nodeTags[i]];
      S = SIN[nodeTags[i]];

      if(_options.debug) {
        #if defined(HAVE_OMP)
        #pragma omp critical
        #endif
        {
          if(_options.smoothDirectionField) {
            fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
            fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -1.5*S, 1.5*C, 0.);
          } else {
            fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
            fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -1.5*S, 1.5*C, 0.);
          }
        }
      }

      double hGrad = 1., hIso = 1.;
      directionsGrad[0] = C;
      directionsGrad[1] = S; 

      getDerivativesAtVertex(_options.polynomialDegree + 1, x, tag, _recoveredFields[0], du, d2u, d3u);

      double normGrad = sqrt(du[0]*du[0] + du[1]*du[1]);
      if(normGrad > 1e-12) {
        computeSizeField(x, tag, directionsGrad, hGrad, hIso, du, d2u, d3u);
      } else {
        // Gradient is zero: region is locally isotropic
        hIso  = _options.hMax;
        hGrad = _options.hMax;
      } 

      double eigenvalues[2] = {fmax(_lambdaMin, fmin(_lambdaMax, 1. / (hGrad * hGrad))), 
                               fmax(_lambdaMin, fmin(_lambdaMax, 1. / (hIso  * hIso )))};
      double ev1[2] = {C, S};
      double ev2[2] = {-S, C};

      MetricTensor M(eigenvalues, ev1, ev2);
      #if defined(HAVE_OMP)
      #pragma omp critical
      #endif
      {
        _metricTensorAtNodetags[nodeTags[i]].assignMatrixFrom(M); // = MetricTensor(eigenvalues, ev1, ev2);
      }

      if(isnan(M(0,0)) || isnan(M(0,1)) || isnan(M(1,0)) || isnan(M(1,1)) ) {
        feInfo("Nan in ISOLINES metric:");
        M.print();
        exit(-1);
      }

      if(M.determinant() < _lambdaMin*_lambdaMin - 1e-10) {
        feInfo("Determinant = %+-1.12e is below lambdaMin*_lambdaMin = %+-1.12e", M.determinant(), _lambdaMin*_lambdaMin);
        exit(-1);
      }
    }

  } // pragma omp parallel

  if(_options.debug) {
    fprintf(directions, "};"); fclose(directions);
    fprintf(dirSmoothed, "};"); fclose(dirSmoothed);
  }

  return FE_STATUS_OK;
}

feStatus feMetric::computeMetricsCurvedLogSimplex(std::vector<std::size_t> &nodeTags,
                                                  std::vector<double> &coord,
                                                  bool useInducedDirections)
{
#if defined(HAVE_SOPLEX)

  size_t numVertices = nodeTags.size();

  // Compute the principal sizes
  const int deg = _options.polynomialDegree;
  const double lMin = _options.hMin;
  const double lMax = _options.hMax;
  const double eps = _options.eTargetError;

  bool OK = true;

  // FILE *directionFile;
  // directionFile = fopen("directions.pos", "w");
  // fprintf(directionFile, "View \" directions \"{\n");

  int cnter = 0;

#if defined(HAVE_OMP)
#pragma omp parallel shared(cnter, OK)
#endif
  {
    double x[2];
    MetricTensor Q;

    int numIter;
    int nTheta = _options.logSimplexOptions.nThetaPerQuadrant;
    int maxIter = _options.logSimplexOptions.maxIter;
    double tol = _options.logSimplexOptions.tol;

    double du[2];

    linearProblem myLP;
    setUpLinearProblem(myLP, _options, nTheta);

// Compute bounded absolute value of upper bound Q at vertices (at nodetags)
#if defined(HAVE_OMP)
#pragma omp for schedule(dynamic, 1)
#endif
    for(size_t i = 0; i < numVertices; i++) {
      if(!OK) continue;

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      getDerivativesAtVertex(1, x, _nodeTag2sequentialTag[nodeTags[i]],  _recoveredFields[0], du);

      // Get the direction of the gradient
      double directionsGrad[2];
      computeDirectionFieldFromGradient(x, _nodeTag2sequentialTag[nodeTags[i]], du, 1e-8, directionsGrad);
      double C = directionsGrad[0];
      double S = directionsGrad[1];
      // #pragma omp critical
      // {
      //   fprintf(directionFile, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0.,  C, S, 0.);
      //   fprintf(directionFile, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.);
      // }

      bool res =
        computeMetricLogSimplexCurved(_nodeTag2sequentialTag[nodeTags[i]], x, directionsGrad,
                                      _recoveredFields[0], Q, maxIter, nTheta, tol, numIter, myLP);

      // If metric computation failed, try again by doubling the number of constraints points in
      // each quadrant (x8 total)
      int retry = 0, maxTry = 5;
      while(!res && retry < maxTry) {
        nTheta *= 2;
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
                   "Could not compute metric. Trying again with %d constraints.", nTheta);
        setUpLinearProblem(myLP, _options, nTheta, true);
        res = computeMetricLogSimplexCurved(_nodeTag2sequentialTag[nodeTags[i]], x, directionsGrad,
                                            _recoveredFields[0], Q, maxIter, nTheta, tol, numIter,
                                            myLP);
        if(res) feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Success");
        retry++;
      }

      if(res) {
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
                   "Computed metric in %2d iterations  - vertex %6d/%6d", numIter, ++cnter,
                   numVertices);

// Bound the eigenvalues of Q
#if defined(HAVE_OMP)
#pragma omp critical
#endif
        {
          if(useInducedDirections) {
            double iso[2] = {-directionsGrad[1], directionsGrad[0]};
            _metricTensorAtNodetags[nodeTags[i]] =
              Q.setEigenvectorsAndBoundEigenvalues(iso, directionsGrad, _lambdaMin, _lambdaMax);
          } else {
            _metricTensorAtNodetags[nodeTags[i]] = Q.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
          }
        }

      } else {
        OK = false;
        // return feErrorMsg(FE_STATUS_ERROR, "Could not compute a metric at (%+-1.5e - %+-1.5e)
        // (vertex %d/%d)",
        //           x[0], x[1], i, nodeTags.size());
      }
    }
  }

  // fprintf(directionFile, "};"); fclose(directionFile);

  if(!OK) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not compute at least one metric tensor");
  }

  return FE_STATUS_OK;
#else
  UNUSED(nodeTags, coord, useInducedDirections);
  return feErrorMsg(
    FE_STATUS_ERROR,
    "SoPlex (optimization library) is required to compute metric tensors for curved P2 adaptation");
#endif
}
void feMetric::averageMetricGradientAtVertices(const std::vector<std::size_t> &nodeTags,
                                               const std::vector<double> &coord,
                                               const std::map<int, MetricTensor> &inputMetrics,
                                               std::map<int, MetricTensor> &dMdx,
                                               std::map<int, MetricTensor> &dMdy)
{
#if defined(HAVE_GMSH)
  // size_t nElm = _recoveredFields[0]->_cnc->getNumElements();
  int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();
  const std::vector<double> &elmArea = _recoveredFields[0]->_cnc->getElementsVolume();

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int elementOrder = getGeometricInterpolantDegree(_recoveredFields[0]->_cnc->getInterpolant());
  int elementType = gmsh::model::mesh::getElementType("Triangle", elementOrder);

  // std::vector<double> weights;
  // gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss5", localCoord, weights);

  // The reference coordinates at which we want the basis functions: just the reference Lagrange nodes here
  std::vector<double> localCoord(3 * nVerticesPerElement, 0.);
  if(nVerticesPerElement == 3) {
    localCoord = {0., 0., 0.,
                  1., 0., 0.,
                  0., 1., 0.};
  } else if(nVerticesPerElement == 6) {
    localCoord = {0.,  0.,  0.,
                  1.,  0.,  0.,
                  0.,  1.,  0.,
                  0.5, 0.,  0.,
                  0.5, 0.5, 0.,
                  0.,  0.5, 0.};
  };

  std::vector<double> basisFunctions;
  std::vector<double> gradBasisFunctions;
  int numComponents, numComponentsGrad, numOrientations, numOrientationsGrad;
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents, basisFunctions, numOrientations);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "GradLagrange", numComponentsGrad, gradBasisFunctions, numOrientationsGrad);
  // int nQuad = weights.size();

  // for(auto val : basisFunctions) {
  //   feInfo("val = %f", val);
  // }

  // Get the mesh elements of the Gmsh model
  std::vector<std::size_t> elementTags;
  std::vector<std::size_t> elemNodeTags;
  gmsh::model::mesh::getElementsByType(elementType, elementTags, elemNodeTags);

  if(elementTags.empty()) {
    feErrorMsg(FE_STATUS_ERROR, "There are no elements in the Gmsh model!");
    exit(-1);
  }

  std::map<int, double> areaAdjacentElements;
  std::map<int, double> m11_x, m12_x, m22_x;
  std::map<int, double> m11_y, m12_y, m22_y;

  size_t numVertices = nodeTags.size();

  for(size_t i = 0; i < numVertices; i++) {

    const MetricTensor &M = inputMetrics.at(nodeTags[i]);
    if(isnan(M(0,0)) || isnan(M(0,1)) || isnan(M(1,0)) || isnan(M(1,1)) ) {
      feInfo("Nan in input metric:");
      M.print();
      exit(-1);
    }

    areaAdjacentElements[nodeTags[i]] = 0.;
    m11_x[nodeTags[i]] = 0.;
    m12_x[nodeTags[i]] = 0.;
    m22_x[nodeTags[i]] = 0.;
    m11_y[nodeTags[i]] = 0.;
    m12_y[nodeTags[i]] = 0.;
    m22_y[nodeTags[i]] = 0.;
  }

  std::vector<double> X(nVerticesPerElement);
  std::vector<double> Y(nVerticesPerElement);

  for(size_t iElm = 0; iElm < elementTags.size(); iElm++) {

    double area = elmArea[iElm];

    int v0 = elemNodeTags[nVerticesPerElement * iElm + 0] - 1;
    int v1 = elemNodeTags[nVerticesPerElement * iElm + 1] - 1;
    int v2 = elemNodeTags[nVerticesPerElement * iElm + 2] - 1;
    int v3, v4, v5;
    if(nVerticesPerElement == 6) {
      v3 = elemNodeTags[nVerticesPerElement * iElm + 3] - 1;
      v4 = elemNodeTags[nVerticesPerElement * iElm + 4] - 1;
      v5 = elemNodeTags[nVerticesPerElement * iElm + 5] - 1;
    }

    double x0 = coord[3 * v0];
    double y0 = coord[3 * v0 + 1];
    double x1 = coord[3 * v1];
    double y1 = coord[3 * v1 + 1];
    double x2 = coord[3 * v2];
    double y2 = coord[3 * v2 + 1];

    double x3, y3, x4, y4, x5, y5;

    if(nVerticesPerElement == 6) {
      x3 = coord[3 * v3]; y3 = coord[3 * v3 + 1];
      x4 = coord[3 * v4]; y4 = coord[3 * v4 + 1];
      x5 = coord[3 * v5]; y5 = coord[3 * v5 + 1];
    }

    X[0] = x0; X[1] = x1; X[2] = x2;
    Y[0] = y0; Y[1] = y1; Y[2] = y2;
    if(nVerticesPerElement == 6) {
      X[3] = x3; X[4] = x4; X[5] = x5;
      Y[3] = y3; Y[4] = y4; Y[5] = y5;
    }

    for(int ii = 0; ii < nVerticesPerElement; ++ii) {

      int ntag_i = elemNodeTags[nVerticesPerElement * iElm + ii];
      areaAdjacentElements.at(ntag_i) += area;

      // Compute inverse of dx/dxsi on element
      // double dxdr[2] = {0., 0.}, dxds[2] = {0., 0.};
      // for(int jj = 0; jj < nVerticesPerElement; ++jj){
      //   double dphi_jdr = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 0];
      //   double dphi_jds = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 1];
      //   dxdr[0] += X[jj] * dphi_jdr;
      //   dxdr[1] += Y[jj] * dphi_jdr;
      //   dxds[0] += X[jj] * dphi_jds;
      //   dxds[1] += Y[jj] * dphi_jds;
      // }

      // double det = dxdr[0]*dxds[1] - dxdr[1]*dxds[0];
      // double drdx[2] = { dxds[1] / det, -dxdr[1] / det};
      // double drdy[2] = {-dxds[0] / det,  dxdr[0] / det};

      double dxdr[2][2] = {{0.,0.},{0.,0.}};
      for(int jj = 0; jj < nVerticesPerElement; ++jj) {
        double dphi_jdr = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 0];
        double dphi_jds = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 1];
        dxdr[0][0] += X[jj] * dphi_jdr;
        dxdr[0][1] += X[jj] * dphi_jds;
        dxdr[1][0] += Y[jj] * dphi_jdr;
        dxdr[1][1] += Y[jj] * dphi_jds;
      }
      double det = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

      double drdx[2][2] = { { dxdr[1][1] / det, -dxdr[0][1] / det}, 
                            {-dxdr[1][0] / det,  dxdr[0][0] / det} };


      double dm11_dr = 0., dm12_dr = 0., dm22_dr = 0.;
      double dm11_ds = 0., dm12_ds = 0., dm22_ds = 0.;

      for(int jj = 0; jj < nVerticesPerElement; ++jj) {
        int ntag_j = elemNodeTags[nVerticesPerElement * iElm + jj];

        double dphi_jdr = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 0];
        double dphi_jds = gradBasisFunctions[3 * nVerticesPerElement * ii + 3 * jj + 1];

        // feInfo("dphidr = %f - %f", dphi_jdr, dphi_jds);

        dm11_dr += inputMetrics.at(ntag_j)(0,0) * dphi_jdr;
        dm12_dr += inputMetrics.at(ntag_j)(0,1) * dphi_jdr;
        dm22_dr += inputMetrics.at(ntag_j)(1,1) * dphi_jdr;

        dm11_ds += inputMetrics.at(ntag_j)(0,0) * dphi_jds;
        dm12_ds += inputMetrics.at(ntag_j)(0,1) * dphi_jds;
        dm22_ds += inputMetrics.at(ntag_j)(1,1) * dphi_jds;
      }

      // m11_x.at(ntag_i) += area * (dm11_dr * drdx[0] + dm11_ds * drdx[1]);
      // m12_x.at(ntag_i) += area * (dm12_dr * drdx[0] + dm12_ds * drdx[1]);
      // m22_x.at(ntag_i) += area * (dm22_dr * drdx[0] + dm22_ds * drdx[1]);

      // m11_y.at(ntag_i) += area * (dm11_dr * drdy[0] + dm11_ds * drdy[1]);
      // m12_y.at(ntag_i) += area * (dm12_dr * drdy[0] + dm12_ds * drdy[1]);
      // m22_y.at(ntag_i) += area * (dm22_dr * drdy[0] + dm22_ds * drdy[1]);

      // feInfo("dMdr = %f - %f - %f - %f - %f - %f", dm11_dr, dm11_ds, dm12_dr, dm12_ds, dm22_dr, dm22_ds);
      // feInfo("drdx = %f - %f - %f - %f", drdx[0][0], drdx[0][1], drdx[1][0], drdx[1][1]);

      m11_x.at(ntag_i) += area * (dm11_dr * drdx[0][0] + dm11_ds * drdx[1][0]);
      m12_x.at(ntag_i) += area * (dm12_dr * drdx[0][0] + dm12_ds * drdx[1][0]);
      m22_x.at(ntag_i) += area * (dm22_dr * drdx[0][0] + dm22_ds * drdx[1][0]);

      m11_y.at(ntag_i) += area * (dm11_dr * drdx[0][1] + dm11_ds * drdx[1][1]);
      m12_y.at(ntag_i) += area * (dm12_dr * drdx[0][1] + dm12_ds * drdx[1][1]);
      m22_y.at(ntag_i) += area * (dm22_dr * drdx[0][1] + dm22_ds * drdx[1][1]);
    }

  }

  for(size_t i = 0; i < numVertices; i++) {
    m11_x.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);
    m12_x.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);
    m22_x.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);
    m11_y.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);
    m12_y.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);
    m22_y.at(nodeTags[i]) /= areaAdjacentElements.at(nodeTags[i]);

    dMdx.at(nodeTags[i])(0,0) = m11_x.at(nodeTags[i]);
    dMdx.at(nodeTags[i])(0,1) = m12_x.at(nodeTags[i]);
    dMdx.at(nodeTags[i])(1,0) = m12_x.at(nodeTags[i]);
    dMdx.at(nodeTags[i])(1,1) = m22_x.at(nodeTags[i]);

    dMdy.at(nodeTags[i])(0,0) = m11_y.at(nodeTags[i]);
    dMdy.at(nodeTags[i])(0,1) = m12_y.at(nodeTags[i]);
    dMdy.at(nodeTags[i])(1,0) = m12_y.at(nodeTags[i]);
    dMdy.at(nodeTags[i])(1,1) = m22_y.at(nodeTags[i]);

    if(isnan(m11_x.at(nodeTags[i])) || isnan(m12_x.at(nodeTags[i])) || isnan(m22_x.at(nodeTags[i])) || isnan(m11_y.at(nodeTags[i])) || isnan(m12_y.at(nodeTags[i])) || isnan(m22_y.at(nodeTags[i])) ) {
      feInfo("Nan:");
      feInfo("Mx = ");
      dMdx.at(nodeTags[i]).print();
      feInfo("My = ");
      dMdy.at(nodeTags[i]).print();
      feInfo("Input metric at node:");
      inputMetrics.at(nodeTags[i]).print();
      feInfo("Determinant = %+-1.10e", inputMetrics.at(nodeTags[i]).determinant());
      feInfo("Adjacent area = %+-1.10e", areaAdjacentElements.at(nodeTags[i]));
      exit(-1);
    }
  }
#else
  UNUSED(nodeTags, coord, inputMetrics, dMdx, dMdy);
#endif
}

//
//
//
//
feStatus feMetric::computeMetricsCurvedReferenceSpace(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (CURVED_REFERENCE_SPACE) for %d vertices...", nodeTags.size());

  size_t numVertices = nodeTags.size();

  // Ugly: create a map between the tags of the patch and GMSH's nodeTags
  // The tags of the patch should be in sequential order.
  double tolNodes = 1e-12;
  std::map<int, int> patchVertices2nodeTags;
  std::vector<int> &patchVertices = _recoveredFields[0]->getVertices();
  for(size_t i = 0; i < patchVertices.size(); i++) {
    int v = patchVertices[i];
    double xv = _recoveredFields[0]->_mesh->getVertex(v)->x();
    double yv = _recoveredFields[0]->_mesh->getVertex(v)->y();
    
    for(size_t j = 0; j < numVertices; j++) {
      const double x = coord[3 * j + 0];
      const double y = coord[3 * j + 1];

      if(fabs(x - xv) < tolNodes && fabs(y - yv) < tolNodes) {
        patchVertices2nodeTags[v] = nodeTags[j];
        break;
      }
    }
  }

  // Create a structure for M^(-1/2) and its derivatives at each vertex
  // The metrics are known at the tags from nodeTags, the GMSH tags.
  std::map<int, MetricTensor> metric_prev;
  std::map<int, MetricTensor> mmud;
  std::map<int, MetricTensor> mmud_dx, new_mmud_dx, new_mmud_dxZN, mmud_dx_prev;
  std::map<int, MetricTensor> mmud_dy, new_mmud_dy, new_mmud_dyZN, mmud_dy_prev;

  MetricTensor Identity(1.0);

  // Initialize M = I, M^(-1/2) = I and its derivatives to zero
  for(size_t i = 0; i < numVertices; i++) {
    _metricTensorAtNodetags[nodeTags[i]].assignMatrixFrom(Identity);
    metric_prev[nodeTags[i]] = MetricTensor(1.0);
    mmud[nodeTags[i]]  = MetricTensor(1.0);

    mmud_dx[nodeTags[i]] = MetricTensor(0.);
    mmud_dy[nodeTags[i]] = MetricTensor(0.);

    mmud_dx_prev[nodeTags[i]] = MetricTensor(0.);
    mmud_dy_prev[nodeTags[i]] = MetricTensor(0.);

    new_mmud_dx[nodeTags[i]] = MetricTensor(0.);
    new_mmud_dy[nodeTags[i]] = MetricTensor(0.);
    new_mmud_dxZN[nodeTags[i]] = MetricTensor(0.);
    new_mmud_dyZN[nodeTags[i]] = MetricTensor(0.);
  }

  #if defined(HAVE_SOPLEX)
  bool applyBinomialCoefficients = true;
  #endif

  bool firstMetric = true;

  bool printMetricConvergence =  _options.curvedMetricReferenceSpaceOptions.printMetricConvergence;

  double damp = _options.curvedMetricReferenceSpaceOptions.dampingCoefficient;
  double alpha = _options.curvedMetricReferenceSpaceOptions.alpha;
  double alphaMax = _options.curvedMetricReferenceSpaceOptions.alphaMax;
  double dalpha = (alphaMax - alpha)/(double) _options.curvedMetricReferenceSpaceOptions.numIterAlpha;

  for(int iAlpha = 0; iAlpha <= _options.curvedMetricReferenceSpaceOptions.numIterAlpha; ++iAlpha, alpha += dalpha) {

    feInfo("Continuation parameter alpha = %+-1.4e", alpha);

    int iterDerivatives = 0;
    bool stopOuter = false;
    bool derivativesConverged = false;

    // Loop until metric derivatives converge
    do {

      // Solve |D^k+1 (f o F(M))(y) zeta^k+1| <= zeta^T * Q * zeta for all zeta with logsimplex method
      //
      // D^k+1 (f o F(M)) is computed with the current metric
      //
      bool OK = true;

      double maxDiffWithPreviousMetric = 0.;

      #if defined(HAVE_OMP)
      #pragma omp parallel shared(OK)
      #endif
      {
        double x[2];
        std::vector<double> errorCoeff(_options.polynomialDegree + 2, 0.);
        double du[2], d2u[4], d3u[8];
        // std::vector<double> DU_EXACT(2, 0.), D2U_EXACT(4, 0.), D3U_EXACT(8, 0.), pos(3, 0.);
        std::vector<double> pos(3, 0.);

        #if defined(HAVE_SOPLEX)
        int numIter;
        int nTheta = _options.logSimplexOptions.nThetaPerQuadrant;
        int maxIter = _options.logSimplexOptions.maxIter;
        double tol = _options.logSimplexOptions.tol;
        linearProblem myLP;
        setUpLinearProblem(myLP, _options, nTheta);
        #endif

        std::vector<double> historicMetricDifferences(_options.curvedMetricReferenceSpaceOptions.maxIterMetrics, 0.);

        // Compute bounded absolute value of upper bound Q at vertices (at nodetags)
        #if defined(HAVE_OMP)
        #pragma omp for schedule(dynamic, 1)
        #endif
        for(size_t i = 0; i < numVertices; i++) {
          if(!OK) continue;

          std::fill(historicMetricDifferences.begin(), historicMetricDifferences.end(), std::nan(""));

          MetricTensor Mprev(1.0), M(1.0);
          M.assignMatrixFrom(_metricTensorAtNodetags.at(nodeTags[i]));
          Mprev.assignMatrixFrom(M);
          MetricTensor Qtri(1.0), Hessian(1.0), Mmud(1.0), Mmud_dx(1.0), Mmud_dy(1.0);

          pos[0] = x[0] = coord[3 * i + 0];
          pos[1] = x[1] = coord[3 * i + 1];
          pos[2] = 0.;

          getDerivativesAtVertex(_options.polynomialDegree + 1, x, _nodeTag2sequentialTag[nodeTags[i]], _recoveredFields[0], du, d2u, d3u);

          // Metric derivatives do not change during computation
          #if defined(HAVE_OMP)
          #pragma omp critical
          #endif
          {
            Mmud_dx.assignMatrixFrom(mmud_dx.at(nodeTags[i]));
            Mmud_dy.assignMatrixFrom(mmud_dy.at(nodeTags[i]));
          }

          if(fabs(Mmud_dx(0,1) - Mmud_dx(1,0)) > 1e-14) { feInfo("Mmud_dx non symmetrique: diff = %+-1.10e", fabs(Mmud_dx(0,1) - Mmud_dx(1,0))); Mmud_dx.print(14); }
          if(fabs(Mmud_dy(0,1) - Mmud_dy(1,0)) > 1e-14) { feInfo("Mmud_dy non symmetrique: diff = %+-1.10e", fabs(Mmud_dy(0,1) - Mmud_dy(1,0))); Mmud_dy.print(14); }

          //
          // Converge the metric at this vertex for fixed metric derivatives
          //
          int iterMetric = 0;
          bool metricIsConverged = false;
          do {

            // SquareMatrix Curvature(2);
            // Curvature(0, 0) = 0.;
            // Curvature(0, 1) = 0.;
            // Curvature(1, 0) = 0.;
            // Curvature(1, 1) = 0.;

            Mmud.assignMatrixFrom(M.pow(-0.5));
            Mmud.symmetrize();
            if(fabs(Mmud(0,1) - Mmud(1,0)) > 1e-14) { feInfo("Mmud non symmetrique: diff = %+-1.10e", fabs(Mmud(0,1) - Mmud(1,0))); Mmud.print(14); exit(-1); }

            if(_options.polynomialDegree == 1) {

              // Error term in the reference space for linear interpolation on quadratic element
              double res[2][2] = {{0., 0.}, {0., 0.}};
              double G[2] = {du[0], du[1]};
              double H[2][2] = {{d2u[0], d2u[1]},
                                {d2u[2], d2u[3]}};

              // For linear interpolation: Compute D^2 (f o F(M))_ij = H_kl Mmud_ki Mmud_lj + G_k dMmud_ki/dx^m Mmud_mj
              for(int I = 0; I < 2; ++I) {
                for(int J = 0; J < 2; ++J) {

                  for(int kk = 0; kk < 2; ++kk) {
                    for(int ll = 0; ll < 2; ++ll) {
                      res[I][J] += H[kk][ll] * Mmud(kk,I) * Mmud(ll,J);
                    }
                  }

                  if(_options.curvedMetricReferenceSpaceOptions.enableCurvature) {
                    for(int kk = 0; kk < 2; ++kk) {
                      for(int mm = 0; mm < 2; ++mm) {
                        if(mm == 0) {
                          // Increment with dMmud/dx
                          res[I][J] += alpha * G[kk] * Mmud_dx(kk,I) * Mmud(mm,J);
                        } else if(mm == 1) {
                          // Increment with dMmud/dy
                          res[I][J] += alpha * G[kk] * Mmud_dy(kk,I) * Mmud(mm,J);
                        } else {
                          feInfo("Error"); exit(-1);
                        }
                      }
                    }
                  }

                  // for(int kk = 0; kk < 2; ++kk) {
                  //   for(int mm = 0; mm < 2; ++mm) {
                  //     if(mm == 0) {
                  //       Curvature(I,J) += alpha * G[kk] * Mmud_dx(kk,I) * Mmud(mm,J);
                  //     } else {
                  //       Curvature(I,J) += alpha * G[kk] * Mmud_dy(kk,I) * Mmud(mm,J);
                  //     }
                  //   }
                  // }

                }
              }

              // // Add a source term for verification
              // if(_options.curvedMetricReferenceSpaceOptions.addSourceTerm) {
              //   if(_options.curvedMetricReferenceSpaceOptions.sourceTerm == nullptr) {
              //     feWarning("Could not add source term because function pointer is null.");
              //     OK = false;
              //   } else {
              //     MetricTensor Source;
              //     _options.curvedMetricReferenceSpaceOptions.sourceTerm->eval(0, pos, Source);
              //     for(int I = 0; I < 2; ++I) {
              //       for(int J = 0; J < 2; ++J) {
              //         res[I][J] += Source(I,J);
              //       }
              //     }
              //   }
              // } 

              if(_options.curvedMetricReferenceSpaceOptions.useAnalyticSolution) {
                // JE NE PENSE PAS QUE ÇA MARCHE PARCE QUE LE TERME COURBE EN PLUS N'EST PAS SYMÉTRIQUE
                // DONC ON DOIT VRAIMENT RÉSOUDRE LE PROBLÈME POUR TROUVER UNE MATRICE SYMÉTRIQUE QUI BORNE
                // LE POLYNOME.

                // SOLUTION: SYMMETRISER LE TERME RES?
                // feInfo("res = %+-1.6e - %+-1.6e - %+-1.6e - %+-1.6e", res[0][0], res[0][1], res[1][0], res[1][1]);
                Hessian(0, 0) = res[0][0];
                Hessian(0, 1) = (res[0][1] + res[1][0]) / 2.;
                Hessian(1, 0) = (res[0][1] + res[1][0]) / 2.;
                Hessian(1, 1) = res[1][1];
              } else {
                errorCoeff[0] =  res[0][0];
                errorCoeff[1] = (res[0][1] + res[1][0]) / 2.;
                errorCoeff[2] =  res[1][1];
              } 

            } else if(_options.polynomialDegree == 2) {

              // Error term in the reference space for quadratic interpolation on quadratic element
              double res[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
              double H[2][2] = {{d2u[0], d2u[1]},
                                {d2u[2], d2u[3]}};
              double C[2][2][2] = {{{d3u[0], d3u[1]}, {d3u[2], d3u[3]}},
                                   {{d3u[4], d3u[5]}, {d3u[6], d3u[7]}}};

              // For quadratic interpolation: Compute D^3 (f o F(M))_ijk = C_mno * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko + ..... LOWER DERIVATIVES
              for(int I = 0; I < 2; ++I) {
                for(int J = 0; J < 2; ++J) {
                  for(int K = 0; K < 2; ++K) {
                    for(int m = 0; m < 2; ++m) {
                      for(int n = 0; n < 2; ++n) {
                        for(int o = 0; o < 2; ++o) {

                          // The term for linear geometry
                          res[I][J][K] += C[m][n][o] * Mmud(I,m) * Mmud(J,n) * Mmud(K,o);

                          if(_options.curvedMetricReferenceSpaceOptions.enableCurvature) {

                            if(o == 0) {
                              res[I][J][K] += alpha * (H[m][n] * Mmud_dx(m,I) * Mmud(o,J) * Mmud(n,K)
                                                     + H[m][n] * Mmud_dx(m,K) * Mmud(o,I) * Mmud(n,J)
                                                     + H[m][n] * Mmud_dx(m,J) * Mmud(o,K) * Mmud(n,I));
                            } else {
                              res[I][J][K] += alpha * (H[m][n] * Mmud_dy(m,I) * Mmud(o,J) * Mmud(n,K)
                                                     + H[m][n] * Mmud_dy(m,K) * Mmud(o,I) * Mmud(n,J)
                                                     + H[m][n] * Mmud_dy(m,J) * Mmud(o,K) * Mmud(n,I));
                            }
                          }

                        }
                      }
                    }
                  }
                }
              }

              if(_options.curvedMetricReferenceSpaceOptions.useAnalyticSolution) {
                errorCoeff[0] =  res[0][0][0];
                errorCoeff[1] =  res[0][0][1] + res[0][1][0] + res[1][0][0];
                errorCoeff[2] =  res[0][1][1] + res[1][0][1] + res[1][1][0];
                errorCoeff[3] =  res[1][1][1];
              } else {
                errorCoeff[0] =  res[0][0][0];
                errorCoeff[1] = (res[0][0][1] + res[0][1][0] + res[1][0][0]) / 3.;
                errorCoeff[2] = (res[0][1][1] + res[1][0][1] + res[1][1][0]) / 3.;
                errorCoeff[3] =  res[1][1][1];
              }

            } else {
              feInfo("TEMPORARY ERROR MESSAGE: CANNOT COMPUTE ANALYTIC DERIVATIVES FOR THIS POLYNOMIAL DEGREE");
              exit(-1);
            }

            /////////////////////////////////////////////////////////////////////////
            // Get reference space upper bound Q_triangle
            bool success = false;
            if(_options.polynomialDegree < 3 && _options.curvedMetricReferenceSpaceOptions.useAnalyticSolution) {

              if(_options.polynomialDegree == 1) {

                Qtri.assignMatrixFrom(Hessian.absoluteValueEigen());
                // Qtri.assignMatrixFrom(Hessian.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax));

                /////////////////////////////////////////
                // Split Hessian and curvature terms
                // MetricTensor foo(1.0);
                // foo(0,0) = Curvature(0,0);
                // foo(0,1) = (Curvature(0,1) + Curvature(1,0)) /2.;
                // foo(1,0) = (Curvature(0,1) + Curvature(1,0)) /2.;
                // foo(1,1) = Curvature(1,1);
                // Qtri.assignMatrixFrom(Hessian.absoluteValueEigen() + foo.absoluteValueEigen());
                /////////////////////////////////////////

                Qtri.symmetrize();
                if(fabs(Qtri(0,1) - Qtri(1,0)) > 1e-14) { feInfo("Qtri non symmetrique: diff = %+-1.10e", fabs(Qtri(0,1) - Qtri(1,0))); Qtri.print(14); exit(-1); }
                success = true;

              } else {
                success = computeAnalyticMetricP2ForLpNorm(errorCoeff, Qtri, 100.);
              } 
            } else {

              #if defined(HAVE_SOPLEX)
              // Solve the linearized minimization problem (i.e., a single iteration of the logsimplex method)
              success = logSimplexOneIteration(x, errorCoeff, _options.polynomialDegree, nTheta, applyBinomialCoefficients, Qtri, myLP);

              // If minimization failed, try again by doubling the number of constraints points in each quadrant (x8 total)
              int retry = 0, maxTry = 5;
              while(!success && retry < maxTry) {
                nTheta *= 2;
                feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Could not compute metric. Trying again with %d constraints.", nTheta);
                setUpLinearProblem(myLP, _options, nTheta, true);
                success = logSimplexOneIteration(x, errorCoeff, _options.polynomialDegree, nTheta, applyBinomialCoefficients, Qtri, myLP);
                if(success) feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Success");
                retry++;
              }
              // Restore initial number of constraints
              nTheta = _options.logSimplexOptions.nThetaPerQuadrant;
              #else
                feWarning("SoPlex (optimization library) is required to compute metric tensors for Pn adaptation");
                OK = false;
              #endif
            }
            /////////////////////////////////////////////////////////////////////////
                   
            /////////////////////////////////////////////////////////////////////////
            // Recover physical space upper bound Q = M^(+1/2) * Q_triangle * M^(+1/2)
            // Update M = Q
            if(success) {

              if(_options.curvedMetricReferenceSpaceOptions.displayQtri) {
                feInfo("Solution for upper bound in reference space: Qtri = %+-1.4e - %+-1.4e - %+-1.4e", Qtri(0,0), Qtri(0,1), Qtri(1,1));
              }

              // Recover Q = M^(1/2) * Qtri * M^(1/2) in the physical space and bound its eigenvalues
              #if defined(HAVE_OMP)
              #pragma omp critical
              #endif
              {
                // MetricTensor Mud = M.sqrt();
                MetricTensor Mud = M.pow(0.5);
                Mud.symmetrize();
                if(fabs(Mud(0,1) - Mud(1,0)) > 1e-14) { feInfo("Mud non symmetrique: diff = %+-1.10e", fabs(Mud(0,1) - Mud(1,0))); Mud.print(14); exit(-1); }
                
                // MetricTensor Q1 = Mud * (Qtri * Mud);
                // MetricTensor Q2 = (Mud * Qtri) * Mud;
                // Q1.symmetrize();
                // Q2.symmetrize();
                // if(fabs(Q1(0,1) - Q1(1,0)) > 1e-14) { feInfo("Q1 non symmetrique: diff = %+-1.10e", fabs(Q1(0,1) - Q1(1,0))); Q1.print(14); }
                // if(fabs(Q2(0,1) - Q2(1,0)) > 1e-14) { feInfo("Q2 non symmetrique: diff = %+-1.10e", fabs(Q2(0,1) - Q2(1,0))); Q2.print(14); }
                // if(fabs(Q1(0,1) - Q1(1,0)) > 1e-14 || fabs(Q2(0,1) - Q2(1,0)) > 1e-14) { exit(-1); }

                MetricTensor Q = (Mud * Qtri) * Mud;
                Q.symmetrize();

                // Update M with Q
                // Damping (useless because we converge the metric independently anyway?)
                if(_options.curvedMetricReferenceSpaceOptions.dampMetrics) {
                  double dampCoeff = _options.curvedMetricReferenceSpaceOptions.dampingMetricCoefficient;
                  M.assignMatrixFrom(M * (1.-dampCoeff) + Q * dampCoeff);
                } else {
                  M.assignMatrixFrom(Q);
                }

                if(fabs(M(0,1) - M(1,0)) > 1e-14) { feInfo("M non symmetrique: diff = %+-1.10e", fabs(M(0,1) - M(1,0))); M.print(14); exit(-1); }

              }
            } else {
              OK = false;
              M(0,0) = _lambdaMin;
              M(0,1) = 0.;
              M(1,0) = 0.;
              M(1,1) = _lambdaMin;
              feWarning("Could not compute a metric at (%+-1.5e - %+-1.5e) (vertex %d/%d)", x[0], x[1],
                         i, nodeTags.size());
            }
            /////////////////////////////////////////////////////////////////////////
            
            /////////////////////////////////////////////////////////////////////////
            // Check convergence and rotate metrics
            double diffMetric = (M - Mprev).frobeniusNorm();
            historicMetricDifferences[iterMetric] = diffMetric;
            if(diffMetric < _options.curvedMetricReferenceSpaceOptions.tolMetrics) {
              metricIsConverged = true;
              if(printMetricConvergence) 
                feInfo("\t\t\t Metric loop at vertex %d - iter %d - diff with previous = %+-1.6e - Converged", i, iterMetric+1, diffMetric);
            } else {
              if(printMetricConvergence)
                feInfo("\t\t\t Metric loop at vertex %d - iter %d - diff with previous = %+-1.6e", i, iterMetric+1, diffMetric);
            }

            // Rotate
            Mprev.assignMatrixFrom(M);
            iterMetric++;
            /////////////////////////////////////////////////////////////////////////
          } while(!metricIsConverged && iterMetric < _options.curvedMetricReferenceSpaceOptions.maxIterMetrics);

          if(!metricIsConverged) {
            OK = false;
            feWarning("Could not compute metric at (%+-1.4e, %+-1.4e) after %d iterations", x[0], x[1], iterMetric);
            feWarning("Historic of last 100 differences with previous metric:");
            for(size_t iHist = historicMetricDifferences.size()-100; iHist < historicMetricDifferences.size(); ++iHist)
              feWarning("Difference history at iter %3d = %+-1.6e", iHist, historicMetricDifferences[iHist]);
          }

          // Metric is converged: assign it in the map
          #if defined(HAVE_OMP)
          #pragma omp critical
          #endif
          {
            // feInfo("Diff with previous M        at vertex %d = %+-1.4e", i, (_metricTensorAtNodetags[nodeTags[i]] - M).frobeniusNorm());
            // feInfo("Diff with previous M^(-1/2) at vertex %d = %+-1.4e", i, (mmud[nodeTags[i]] - M.pow(-0.5)).frobeniusNorm());
            double diff = (_metricTensorAtNodetags.at(nodeTags[i]) - M).frobeniusNorm();
            maxDiffWithPreviousMetric = fmax(maxDiffWithPreviousMetric, diff);

            _metricTensorAtNodetags.at(nodeTags[i]).assignMatrixFrom(M);
            mmud.at(nodeTags[i]).assignMatrixFrom(M.pow(-0.5));
          }

        } // for numVertices

      } // omp parallel

      if(!OK) {
        return feErrorMsg(FE_STATUS_ERROR, "Could not compute at least one metric tensor");
      }

      // feInfo("\t\tMetric converged. Max difference with previous metric = %+-1.4e", maxDiffWithPreviousMetric);
      // Plot always first metric field
      if(firstMetric) {
        firstMetric = false;
        drawEllipsoids("metrics_first.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 100);
      }

      if(_options.curvedMetricReferenceSpaceOptions.debugPlot) {
        drawEllipsoids("metrics_iAlpha_" + std::to_string(iAlpha) + "_iter" + std::to_string(iterDerivatives) + ".pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
      }

      if(_options.curvedMetricReferenceSpaceOptions.enableGradation) {
        newGradation(nodeTags, coord, _options.curvedMetricReferenceSpaceOptions.gradation, _metricTensorAtNodetags);
        if(_options.curvedMetricReferenceSpaceOptions.debugPlot) {
          drawEllipsoids("metrics_iAlpha_" + std::to_string(iAlpha) + "_iter" + std::to_string(iterDerivatives) + "_afterGradation.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        }
      }

      // Recompute derivatives of M^(-1/2)
      if(_options.curvedMetricReferenceSpaceOptions.averageMetricGradient) {
        averageMetricGradientAtVertices(nodeTags, coord, mmud, new_mmud_dx, new_mmud_dy);
        // _recoveredFields[0]->computeRHSAndSolve_noIntegral_inputMetric(patchVertices2nodeTags, mmud, new_mmud_dxZN, new_mmud_dyZN);
      } else {
        // Recover the gradient of each component with Zhang and Naga's method
        _recoveredFields[0]->computeRHSAndSolve_noIntegral_inputMetric(patchVertices2nodeTags, mmud, new_mmud_dx, new_mmud_dy);
      }
      
      double maxDiffDerivatives = 0.;
      for(size_t i = 0; i < numVertices; i++) {

        // feInfo(" ");
        // feInfo("dMmud_dx from averaging = %+-1.6e - %+-1.6e - %+-1.6e - %+-1.6e",
        //   new_mmud_dx.at(nodeTags[i])(0,0), new_mmud_dx.at(nodeTags[i])(0,1), new_mmud_dx.at(nodeTags[i])(1,0), new_mmud_dx.at(nodeTags[i])(1,1));
        // feInfo("dMmud_dx from ZN        = %+-1.6e - %+-1.6e - %+-1.6e - %+-1.6e",
        //   new_mmud_dxZN.at(nodeTags[i])(0,0), new_mmud_dxZN.at(nodeTags[i])(0,1), new_mmud_dxZN.at(nodeTags[i])(1,0), new_mmud_dxZN.at(nodeTags[i])(1,1));
        // feInfo("dMmud_dy from averaging = %+-1.6e - %+-1.6e - %+-1.6e - %+-1.6e",
        //   new_mmud_dy.at(nodeTags[i])(0,0), new_mmud_dy.at(nodeTags[i])(0,1), new_mmud_dy.at(nodeTags[i])(1,0), new_mmud_dy.at(nodeTags[i])(1,1));
        // feInfo("dMmud_dy from ZN        = %+-1.6e - %+-1.6e - %+-1.6e - %+-1.6e",
        //   new_mmud_dyZN.at(nodeTags[i])(0,0), new_mmud_dyZN.at(nodeTags[i])(0,1), new_mmud_dyZN.at(nodeTags[i])(1,0), new_mmud_dyZN.at(nodeTags[i])(1,1));

        // Damping
        if(_options.curvedMetricReferenceSpaceOptions.dampDerivatives) {
          mmud_dx.at(nodeTags[i]).assignMatrixFrom(mmud_dx_prev.at(nodeTags[i]) * (1.-damp) + new_mmud_dx.at(nodeTags[i]) * damp);
          mmud_dy.at(nodeTags[i]).assignMatrixFrom(mmud_dy_prev.at(nodeTags[i]) * (1.-damp) + new_mmud_dy.at(nodeTags[i]) * damp);
        } else {
          mmud_dx.at(nodeTags[i]).assignMatrixFrom(new_mmud_dx.at(nodeTags[i]));
          mmud_dy.at(nodeTags[i]).assignMatrixFrom(new_mmud_dy.at(nodeTags[i]));
        }

        // Check convergence for outer loop
        double diffx = (mmud_dx_prev.at(nodeTags[i]) - mmud_dx.at(nodeTags[i])).frobeniusNorm();
        double diffy = (mmud_dy_prev.at(nodeTags[i]) - mmud_dy.at(nodeTags[i])).frobeniusNorm();
        maxDiffDerivatives = fmax(maxDiffDerivatives, diffx);
        maxDiffDerivatives = fmax(maxDiffDerivatives, diffy);

        if(_options.curvedMetricReferenceSpaceOptions.displayDerivatives) {
          feInfo("Computed dx = %+-1.4e - %+-1.4e - %+-1.4e - Diff = %+-1.4e", mmud_dx[nodeTags[i]](0,0), mmud_dx[nodeTags[i]](0,1), mmud_dx[nodeTags[i]](1,1), diffx);
          feInfo("Computed dy = %+-1.4e - %+-1.4e - %+-1.4e - Diff = %+-1.4e", mmud_dy[nodeTags[i]](0,0), mmud_dy[nodeTags[i]](0,1), mmud_dy[nodeTags[i]](1,1), diffy);
          feInfo(" ");
        }
      }
      // if(iterDerivatives > 1)
      //   exit(-1);

      if(maxDiffDerivatives < _options.curvedMetricReferenceSpaceOptions.tolDerivatives) {
        derivativesConverged = true;
        // feInfo("\t\tExiting derivatives loop because difference w.r.t. previous derivatives is below tolerance: %+-1.6e < %+-1.6e",
        //   maxDiffDerivatives, _options.curvedMetricReferenceSpaceOptions.tolDerivatives);
        feInfo("\t\tIter %2d - All metrics converged with max change = %+-1.4e - Current change on derivatives = %+-1.4e - Converged for tol = %+-1.4e",
          iterDerivatives, maxDiffWithPreviousMetric, maxDiffDerivatives,  _options.curvedMetricReferenceSpaceOptions.tolDerivatives);
      } else{
        feInfo("\t\tIter %2d - All metrics converged with max change = %+-1.4e - Current change on derivatives = %+-1.4e",
          iterDerivatives, maxDiffWithPreviousMetric, maxDiffDerivatives);
      }

      // Rotate the metrics derivatives
      for(size_t i = 0; i < numVertices; i++) {
        mmud_dx_prev.at(nodeTags[i]).assignMatrixFrom(mmud_dx.at(nodeTags[i]));
        mmud_dy_prev.at(nodeTags[i]).assignMatrixFrom(mmud_dy.at(nodeTags[i]));
      }

      iterDerivatives++;

      stopOuter = false;

      // if(iterDerivatives <= _options.curvedMetricReferenceSpaceOptions.numIterAlpha) {
      //   // Increase alpha
      //   alpha += dalpha;
      // } else {
        if(derivativesConverged) {
          stopOuter = true;
        }
      // }

      if(iterDerivatives > _options.curvedMetricReferenceSpaceOptions.maxIterDerivatives) {
        if(iterDerivatives <= _options.curvedMetricReferenceSpaceOptions.numIterAlpha) {
          feWarning("Exiting outer loop for derivatives but continuation parameter did not reach 1.");
          feWarning("Last converged value of continuation parameter = %+-1.4e", alpha - dalpha);
        }
        stopOuter = true;
      }

    } while(!stopOuter);

    // Plot always last metric field
    drawEllipsoids("metrics_last.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 100);

    if(!derivativesConverged) {
      return feErrorMsg(FE_STATUS_ERROR, "Metric derivatives did not converge after %d iterations", iterDerivatives);
    }

  } // for iAlpha

  return FE_STATUS_OK;
}

//
// Compute the error commited when computing the derivatives of M^(-1/2)
// For development and debug only
// The metric derivatives are interpolated using simply the derivatives of the shape functions
// and are compared with exact derivatives or with finite difference derivatives
// 
void feMetric::computeErrorOnMetricDerivatives(double &errorOnMmud_dx, double &errorOnMmud_dy)
{
#if defined(HAVE_GMSH)

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  feStatus s = createVertex2NodeMap(nodeTags, coord);
  if(s != FE_STATUS_OK) {
    feWarning("Could not create Gmsh node tags");
    exit(-1);
  }

  int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  if(elementTypes.size() > 1) {
    feErrorMsg(FE_STATUS_ERROR, "Gmsh mesh has more than one element type."
                                "We expect only P1 or P2 triangles in 2D entities");
    exit(-1);
  }

  int elementType = elementTypes[0];

  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss6", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);
  std::vector<double> gradBasisFunctions;
  int numComponentsGrad, numOrientationsGrad;
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "GradLagrange", numComponentsGrad, gradBasisFunctions, numOrientationsGrad);

  if(_nVerticesPerElmOnBackmesh == 3 && nVerticesPerElement == 6) {
    // Linear backmesh with a P2 feMesh : get P1 basis
    int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
    gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents,
                                         basisFunctions, numOrientations);
  }

  // Get the jacobians
  std::vector<double> jac, detJac, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, detJac, pts);

  double I = 0.0;
  // double E_DROIT = 0.0;
  // double E_CURVE = 0.0;
  // double COMP = 0.0;

  double hFD = 1e-6;

  size_t nQuad = weights.size();
  size_t nElm = elementTags[0].size();

  size_t numVertices = nodeTags.size();

  // Compute M^(-1/2) and its derivatives
  for(size_t i = 0; i < numVertices; i++) {
    _mmud[nodeTags[i]].assignMatrixFrom(_metricTensorAtNodetags.at(nodeTags[i]).pow(-0.5));
    _dmmud_dx[nodeTags[i]] = MetricTensor(0.0);
    _dmmud_dy[nodeTags[i]] = MetricTensor(0.0);
  }

  // Ugly: create a map between the tags of the patch and GMSH's nodeTags
  // The tags of the patch should be in sequential order.
  double tol = 1e-12;
  std::map<int, int> patchVertices2nodeTags;
  std::vector<int> &patchVertices = _recoveredFields[0]->getVertices();
  for(size_t i = 0; i < patchVertices.size(); i++) {
    int v = patchVertices[i];
    double xv = _recoveredFields[0]->_mesh->getVertex(v)->x();
    double yv = _recoveredFields[0]->_mesh->getVertex(v)->y();
    
    for(size_t j = 0; j < numVertices; j++) {
      const double x = coord[3 * j + 0];
      const double y = coord[3 * j + 1];

      if(fabs(x - xv) < tol && fabs(y - yv) < tol) {
        patchVertices2nodeTags[v] = nodeTags[j];
        break;
      }
    }
  }

  // Compute derivatives with Zhang/Naga or by averaging the FE gradients at vertices
  // _recoveredFields[0]->computeRHSAndSolve_noIntegral_inputMetric(patchVertices2nodeTags, _mmud, _dmmud_dx, _dmmud_dy);
  averageMetricGradientAtVertices(nodeTags, coord, _mmud, _dmmud_dx, _dmmud_dy);

  double E_DMMUD_DX = 0.;
  double E_DMMUD_DY = 0.;

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double xsi[2];
    feFunctionArguments args;
    MetricTensor M, fooM, M_interp, Mmud_ref, dMdx_ref, dMdx, dMmud_dx, dMmud_dy, dMmud_dx_ref, dMmud_dy_ref;
    std::vector<MetricTensor> MX_interpolation(_nVerticesPerElmOnBackmesh);
    std::vector<MetricTensor> MY_interpolation(_nVerticesPerElmOnBackmesh);
    std::vector<MetricTensor> Mmud_at_nodes(_nVerticesPerElmOnBackmesh);
    std::vector<MetricTensor> M_at_nodes(_nVerticesPerElmOnBackmesh);
    std::vector<MetricTensor> logM_at_nodes(_nVerticesPerElmOnBackmesh);

    feFunctionArguments argsFD;
    MetricTensor mfd, mmud_ref_FD, mx_ref_FD, my_ref_FD, mfdx, mfdy;

    /////////////////////////////////////////////////////////
    std::vector<double> X(_nVerticesPerElmOnBackmesh);
    std::vector<double> Y(_nVerticesPerElmOnBackmesh);
    /////////////////////////////////////////////////////////

#if defined(HAVE_OMP)
#pragma omp for reduction(+ : I, E_DMMUD_DX, E_DMMUD_DY)
#endif
    for(size_t iElm = 0; iElm < nElm; iElm++) {

      /////////////////////////////////////////////////////////
      // Get vertices coordinates
      int v0 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 0] - 1;
      int v1 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 1] - 1;
      int v2 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 2] - 1;
      int v3, v4, v5;
      if(_nVerticesPerElmOnBackmesh == 6) {
        v3 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 3] - 1;
        v4 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 4] - 1;
        v5 = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 5] - 1;
      }

      double x0 = coord[3 * v0];
      double y0 = coord[3 * v0 + 1];
      double x1 = coord[3 * v1];
      double y1 = coord[3 * v1 + 1];
      double x2 = coord[3 * v2];
      double y2 = coord[3 * v2 + 1];
      double x3, y3, x4, y4, x5, y5;
      if(_nVerticesPerElmOnBackmesh == 6) {
        x3 = coord[3 * v3]; y3 = coord[3 * v3 + 1];
        x4 = coord[3 * v4]; y4 = coord[3 * v4 + 1];
        x5 = coord[3 * v5]; y5 = coord[3 * v5 + 1];
      }
      X[0] = x0; X[1] = x1; X[2] = x2;
      Y[0] = y0; Y[1] = y1; Y[2] = y2;
      if(_nVerticesPerElmOnBackmesh == 6) {
        X[3] = x3; X[4] = x4; X[5] = x5;
        Y[3] = y3; Y[4] = y4; Y[5] = y5;
      }
      /////////////////////////////////////////////////////////


      for(int i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
        int tag = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + i];
        MX_interpolation[i].assignMatrixFrom(_dmmud_dx.at(tag));
        MY_interpolation[i].assignMatrixFrom(_dmmud_dy.at(tag));

        std::vector<double> pos{coord[3 * (tag-1) + 0], coord[3 * (tag-1) + 1], 0.};
        args.pos[0] = pos[0];
        args.pos[1] = pos[1];
        args.pos[2] = pos[2];
        _options.analyticMetric->eval(args, fooM);  
        Mmud_at_nodes[i].assignMatrixFrom(fooM.pow(-0.5));
        M_at_nodes[i].assignMatrixFrom(fooM);
        logM_at_nodes[i].assignMatrixFrom(fooM.log());
      }

      for(size_t k = 0; k < nQuad; k++) {

        xsi[0] = localCoord[3 * k + 0];
        xsi[1] = localCoord[3 * k + 1];

        // Coordonnees physiques
        args.pos[0] = pts[iElm * 3 * nQuad + 3 * k + 0];
        args.pos[1] = pts[iElm * 3 * nQuad + 3 * k + 1];
        args.pos[2] = 0.;
        
        // Evaluate analytic M^(-1/2) and derivatives
        _options.analyticMetric->eval(args, M);
        Mmud_ref.assignMatrixFrom(M.pow(-0.5));
        _options.analytic_dMmud_dx->eval(args, dMmud_dx_ref);
        _options.analytic_dMmud_dy->eval(args, dMmud_dy_ref);

        // // Interpolate derivatives of M^(-1/2) at quadrature node
        // // Derivatives are not SPD, so cant use log-euclidian interpolation
        // if(_nVerticesPerElmOnBackmesh == 3) {
        //   dMmud_dx.assignMatrixFrom(MX_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
        //                           + MX_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
        //                           + MX_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]);
        //   dMmud_dy.assignMatrixFrom(MY_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
        //                           + MY_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
        //                           + MY_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]);
        // } else if(_nVerticesPerElmOnBackmesh == 6) {
        //   dMmud_dx.assignMatrixFrom(MX_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
        //                           + MX_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
        //                           + MX_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]
        //                           + MX_interpolation[3] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 3]
        //                           + MX_interpolation[4] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 4]
        //                           + MX_interpolation[5] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 5]);
        //   dMmud_dy.assignMatrixFrom(MY_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
        //                           + MY_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
        //                           + MY_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]
        //                           + MY_interpolation[3] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 3]
        //                           + MY_interpolation[4] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 4]
        //                           + MY_interpolation[5] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 5]);
        // } else {
        //   feInfo("Number of vertices per element on backmesh not supported. Expected P1 or P2 triangle.");
        //   exit(-1);
        // }

        // OR
        // interpolate M^(-1/2) with the gradient of the basis functions
        if(_nVerticesPerElmOnBackmesh == 3) {

          /////////////////////////////////////////////////////////
          // Compute inverse of dx/dxsi on element
          double dxdr[2][2] = {{0.,0.},{0.,0.}};
          for(int jj = 0; jj < nVerticesPerElement; ++jj) {
            double dphi_jdr = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * jj + 0];
            double dphi_jds = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * jj + 1];
            dxdr[0][0] += X[jj] * dphi_jdr;
            dxdr[0][1] += X[jj] * dphi_jds;
            dxdr[1][0] += Y[jj] * dphi_jdr;
            dxdr[1][1] += Y[jj] * dphi_jds;
          }
          double det = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];

          double drdx[2][2] = { { dxdr[1][1] / det, -dxdr[0][1] / det}, 
                                {-dxdr[1][0] / det,  dxdr[0][0] / det} };

          // double phi0 = basisFunctions[_nVerticesPerElmOnBackmesh * k + 0];
          // double phi1 = basisFunctions[_nVerticesPerElmOnBackmesh * k + 1];
          // double phi2 = basisFunctions[_nVerticesPerElmOnBackmesh * k + 2];

          double dphi0_jdr = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 0 + 0];
          double dphi1_jdr = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 1 + 0];
          double dphi2_jdr = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 2 + 0];
          double dphi0_jds = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 0 + 1];
          double dphi1_jds = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 1 + 1];
          double dphi2_jds = gradBasisFunctions[3 * _nVerticesPerElmOnBackmesh * k + 3 * 2 + 1];

          // double sumPhir = dphi0_jdr + dphi1_jdr + dphi2_jdr;
          // double sumPhis = dphi0_jds + dphi1_jds + dphi2_jds;
          // feInfo("sumPhir = %f - r = %f", sumPhir, sumPhis);

          double dphi0_dx = dphi0_jdr * drdx[0][0] + dphi0_jds * drdx[1][0];
          double dphi1_dx = dphi1_jdr * drdx[0][0] + dphi1_jds * drdx[1][0];
          double dphi2_dx = dphi2_jdr * drdx[0][0] + dphi2_jds * drdx[1][0];

          /////////////////////////////////////////////////////////
          // M_interp.assignMatrixFrom(M_at_nodes[0] * phi0 + M_at_nodes[1] * phi1 + M_at_nodes[2] * phi2);
          // logEuclidianP1Interpolation(xsi, logM_at_nodes[0], logM_at_nodes[1], logM_at_nodes[2], M_interp);

          // Derivatives of M
          dMdx.assignMatrixFrom(M_at_nodes[0] * dphi0_dx
                              + M_at_nodes[1] * dphi1_dx
                              + M_at_nodes[2] * dphi2_dx);

          // Derivatives of M^(-1/2) using gradient of basis functions
          // dMmud_dx.assignMatrixFrom(Mmud_at_nodes[0] * dphi0_dx
          //                         + Mmud_at_nodes[1] * dphi1_dx
          //                         + Mmud_at_nodes[2] * dphi2_dx);

          // Derivatives of M^(-1/2) using dM^(-1/2) at nodes and the basis functions
          dMmud_dx.assignMatrixFrom(MX_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
                                  + MX_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
                                  + MX_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]);

          dMmud_dy.assignMatrixFrom(MY_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
                                  + MY_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
                                  + MY_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * k + 2]);

        } else if(_nVerticesPerElmOnBackmesh == 6) {

          logEuclidianP2Interpolation(xsi, logM_at_nodes[0], logM_at_nodes[1], logM_at_nodes[2], 
                                           logM_at_nodes[3], logM_at_nodes[4], logM_at_nodes[5], M_interp);

          // FAUX: changer pour le gradient physique ou interpoler les dérivées aux noeuds comme ci-dessus
          // dMmud_dx.assignMatrixFrom(Mmud_at_nodes[0] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
          //                         + Mmud_at_nodes[1] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
          //                         + Mmud_at_nodes[2] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 2]
          //                         + Mmud_at_nodes[3] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 3]
          //                         + Mmud_at_nodes[4] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 4]
          //                         + Mmud_at_nodes[5] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 5]);
          // dMmud_dy.assignMatrixFrom(Mmud_at_nodes[0] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 0]
          //                         + Mmud_at_nodes[1] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 1]
          //                         + Mmud_at_nodes[2] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 2]
          //                         + Mmud_at_nodes[3] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 3]
          //                         + Mmud_at_nodes[4] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 4]
          //                         + Mmud_at_nodes[5] * gradBasisFunctions[_nVerticesPerElmOnBackmesh * k + 5]);
        } else {
          feInfo("Number of vertices per element on backmesh not supported. Expected P1 or P2 triangle.");
          exit(-1);
        }

        // Interpolation of M
        // double diff_M = (M - M_interp).frobeniusNorm();
        // double diff_M = matNorm2(M, M_interp);
        // feInfo("M        = %+-1.4e - %+-1.4e - %+-1.4e", M(0,0), M(0,1), M(1,1));
        // feInfo("M_interp = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e", M_interp(0,0), M_interp(0,1), M_interp(1,1), diff_M);
        // I += weights[k] * detJac[iElm * nQuad + k];
        // E_DMMUD_DX += weights[k] * detJac[iElm * nQuad + k] * (diff_M*diff_M);
        // E_DMMUD_DY += weights[k] * detJac[iElm * nQuad + k] * (0);

        // // Interpolation of dMdx (hardcoded analytic solution): OK
        // double a = 1.234567;
        // double b = 0.012345;
        // double c = 5.987654;
        // dMdx_ref(0,0) = 2.*a*pos[0];
        // dMdx_ref(0,1) = 2.*b*pos[0];
        // dMdx_ref(1,0) = 2.*b*pos[0];
        // dMdx_ref(1,1) = 0.;
        // double diff_dMdx = matNorm2(dMdx_ref, dMdx);
        // feInfo("dMdx        = %+-1.4e - %+-1.4e - %+-1.4e", dMdx_ref(0,0), dMdx_ref(0,1), dMdx_ref(1,1));
        // feInfo("dMdx_interp = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e", dMdx(0,0), dMdx(0,1), dMdx(1,1), diff_dMdx);
        // I += weights[k] * detJac[iElm * nQuad + k];
        // E_DMMUD_DX += weights[k] * detJac[iElm * nQuad + k] * (diff_dMdx*diff_dMdx);
        // E_DMMUD_DY += weights[k] * detJac[iElm * nQuad + k] * (0);

        // Interpolation of dM^(-1/2)/dx
        double diff_Mx = matNorm2(dMmud_dx_ref, dMmud_dx);
        double diff_My = matNorm2(dMmud_dy_ref, dMmud_dy);
        feInfo("Mmudx        = %+-1.4e - %+-1.4e - %+-1.4e", dMmud_dx_ref(0,0), dMmud_dx_ref(0,1), dMmud_dx_ref(1,1));
        feInfo("Mmudx_interp = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e", dMmud_dx(0,0), dMmud_dx(0,1), dMmud_dx(1,1), diff_Mx);
        argsFD.pos[0] = args.pos[0] + hFD;
        argsFD.pos[1] = args.pos[1];
        argsFD.pos[2] = args.pos[2];
        _options.analyticMetric->eval(argsFD, mfd);
        mmud_ref_FD.assignMatrixFrom(mfd.pow(-0.5));
        mfdx.assignMatrixFrom((mmud_ref_FD - Mmud_ref) * (1./hFD));        
        double diff_FD = matNorm2(dMmud_dx_ref, mfdx);
        double diff_FD_interp = matNorm2(dMmud_dx, mfdx);
        feInfo("         FD  = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - diff with interp = %+-1.4e", mfdx(0,0), mfdx(0,1), mfdx(1,1), diff_FD, diff_FD_interp);

        I += weights[k] * detJac[iElm * nQuad + k];
        E_DMMUD_DX += weights[k] * detJac[iElm * nQuad + k] * (diff_Mx*diff_Mx);
        E_DMMUD_DY += weights[k] * detJac[iElm * nQuad + k] * (diff_My*diff_My);

        // double diff_dx = (dMmud_dx - dMmud_dx_ref).frobeniusNorm();
        // double diff_dy = (dMmud_dy - dMmud_dy_ref).frobeniusNorm();
        
        // ///////////////////////////////////////////////////////////////
        // // Differences finies
        // posFD[0] = pos[0] + hFD;
        // posFD[1] = pos[1];
        // posFD[2] = pos[2];
        // _options.analyticMetric->eval(argsFD, mfd);
        // mmud_ref_FD.assignMatrixFrom(mfd.pow(-0.5));
        // mfdx.assignMatrixFrom((mmud_ref_FD - Mmud_ref) * (1./hFD));
        // double diff_FD_x = (mfdx - dMmud_dx_ref).frobeniusNorm();

        // feInfo("dx ref  = %+-1.4e - %+-1.4e - %+-1.4e", dMmud_dx_ref(0,0), dMmud_dx_ref(0,1), dMmud_dx_ref(1,1));
        // feInfo(" interp = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e", dMmud_dx(0,0), dMmud_dx(0,1), dMmud_dx(1,1), diff_dx);
        // feInfo("    FD  = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e", mfdx(0,0), mfdx(0,1), mfdx(1,1), diff_FD_x);

        // ///////////////////////////////////////////////////////////////


        // I += weights[k] * detJac[iElm * nQuad + k];
        // E_DMMUD_DX += weights[k] * detJac[iElm * nQuad + k] * (diff_dx*diff_dx);
        // E_DMMUD_DY += weights[k] * detJac[iElm * nQuad + k] * (diff_dy*diff_dy);

      } // for quad nodes
    } // for elm
  } // omp parallel

  // // ESSAI: erreur L2 discrete sur les derivees aux sommets
  // std::vector<double> mypos(3,0.);
  // MetricTensor foommudx, nodal;
  // E_DMMUD_DX = 0.;
  // for(size_t i = 0; i < numVertices; ++i) {
  //   int tag = nodeTags[i];
  //   mypos[0] = coord[3 * (tag-1) + 0];
  //   mypos[1] = coord[3 * (tag-1) + 1];
  //   _options.analytic_dMmud_dx->eval(0, mypos, foommudx);
  //   MetricTensor &nodal = _dmmud_dx.at(tag);
  //   double err = matNorm2(nodal, foommudx);
  //   double err_rel = err/foommudx.frobeniusNorm();
  //   feInfo("Mmudx        = %+-1.4e - %+-1.4e - %+-1.4e", foommudx(0,0), foommudx(0,1), foommudx(1,1));
  //   feInfo("Mmudx nodal  = %+-1.4e - %+-1.4e - %+-1.4e - diff = %+-1.4e - rel = %1.4e", nodal(0,0), nodal(0,1), nodal(1,1), err, err_rel);
  //   E_DMMUD_DX += err*err;

  //   if(foommudx.frobeniusNorm() > 1e-10 && err_rel > 10.) {
  //     feInfo("Big error at %f - %f", mypos[0], mypos[1]);
  //     exit(-1);
  //   }
  // }

  E_DMMUD_DX = sqrt(E_DMMUD_DX);
  E_DMMUD_DY = sqrt(E_DMMUD_DY);

  feInfoCond(FE_VERBOSE > 0, "Integral = %+-1.16e", I);
  feInfoCond(FE_VERBOSE > 0, "E_DMMUD_DX = %+-1.16e", E_DMMUD_DX);
  feInfoCond(FE_VERBOSE > 0, "E_DMMUD_DY = %+-1.16e", E_DMMUD_DY);

  errorOnMmud_dx = E_DMMUD_DX;
  errorOnMmud_dy = E_DMMUD_DY;
#else
  UNUSED(errorOnMmud_dx, errorOnMmud_dy);
#endif
}

//
// Compute the continuous error model for this metric and the precribed field,
// which can be given through analytic or reconstructed derivatives.
// 
void feMetric::computeContinuousErrorModel(const bool withAnalyticMetric, bool computeError,
  bool plotError, const std::string &errorFileName, double &errorStraight, double &errorCurved, double &complexity)
{
#if defined(HAVE_GMSH)

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  feStatus s = createVertex2NodeMap(nodeTags, coord);
  if(s != FE_STATUS_OK) {
    feWarning("Could not create Gmsh node tags");
    exit(-1);
  }

  int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  if(elementTypes.size() > 1) {
    feErrorMsg(FE_STATUS_ERROR, "Gmsh mesh has more than one element type."
                                "We expect only P1 or P2 triangles in 2D entities");
    exit(-1);
  }

  int elementType = elementTypes[0];

  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss20", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents,
                                       basisFunctions, numOrientations);

  if(_nVerticesPerElmOnBackmesh == 3 && nVerticesPerElement == 6) {
    // Linear backmesh with a P2 feMesh : get P1 basis
    int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
    gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents,
                                         basisFunctions, numOrientations);
  }

  // Get the jacobians
  std::vector<double> jac, detJac, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, detJac, pts);

  size_t numVertices = nodeTags.size();

  // Compute M^(-1/2) and its derivatives
  for(size_t i = 0; i < numVertices; i++) {
    MetricTensor &M = _metricTensorAtNodetags.at(nodeTags[i]);
    MetricTensor myMmud = M.pow(-0.5);
    if(isnan(myMmud(0,0)) || isnan(myMmud(0,1)) || isnan(myMmud(1,0)) || isnan(myMmud(1,1)) ) {
      feInfo("Nan in inverse square root metric:");
      M.print(14);
      feInfo("det M = %+-1.10e", M.determinant());
      myMmud.print(14);
      exit(-1);
    }
    _mmud[nodeTags[i]].assignMatrixFrom(myMmud);
    _dmmud_dx[nodeTags[i]] = MetricTensor(0.0);
    _dmmud_dy[nodeTags[i]] = MetricTensor(0.0);
  }

  if(!withAnalyticMetric) {

    // Ugly: create a map between the tags of the patch and GMSH's nodeTags
    // The tags of the patch should be in sequential order.
    double tol = 1e-12;
    std::map<int, int> patchVertices2nodeTags;
    std::vector<int> &patchVertices = _recoveredFields[0]->getVertices();
    for(size_t i = 0; i < patchVertices.size(); i++) {
      int v = patchVertices[i];
      double xv = _recoveredFields[0]->_mesh->getVertex(v)->x();
      double yv = _recoveredFields[0]->_mesh->getVertex(v)->y();
      
      for(size_t j = 0; j < numVertices; j++) {
        const double x = coord[3 * j + 0];
        const double y = coord[3 * j + 1];

        if(fabs(x - xv) < tol && fabs(y - yv) < tol) {
          patchVertices2nodeTags[v] = nodeTags[j];
          break;
        }
      }
    }

    // Compute derivatives with Zhang/Naga or by averaging the FE gradients at vertices
    // _recoveredFields[0]->computeRHSAndSolve_noIntegral_inputMetric(patchVertices2nodeTags, _mmud, _dmmud_dx, _dmmud_dy); // ZN
    averageMetricGradientAtVertices(nodeTags, coord, _mmud, _dmmud_dx, _dmmud_dy); // Averaging
  }

  double Area = 0.0;
  double E_DROIT = 0.0;
  double E_CURVE = 0.0;
  double COMP = 0.0;

  size_t nQuad = weights.size();
  size_t nElm = elementTags[0].size();

  std::string linearErrorFile = errorFileName + "_linear.pos";
  std::string curvedErrorFile = errorFileName + "_curved.pos";

  FILE *fileLinear = nullptr, *fileCurved = nullptr;
  if(plotError) {
    fileLinear = fopen(linearErrorFile.data(), "w");
    fileCurved = fopen(curvedErrorFile.data(), "w");
    fprintf(fileLinear, "View\" %s \"{\n", linearErrorFile.data());
    fprintf(fileCurved, "View\" %s \"{\n", curvedErrorFile.data());
  }

  std::vector<double> errorLinearOnElements(nElm, 0.);
  std::vector<double> errorCurvedOnElements(nElm, 0.);

#define MAX_VERTICES_BACKMESH 6

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    double G[2];
    MetricTensor absH, H, Eloc, Eloc_curved, M, dMmud_dx, dMmud_dy, Mmud;
    MetricTensor logM0, logM1, logM2, logM3, logM4, logM5;
    MetricTensor Qtri;
    double xsi[2];
    int tags[MAX_VERTICES_BACKMESH];
    std::vector<double> DU_EXACT(2, 0.), D2U_EXACT(4, 0.), D3U_EXACT(8, 0.);
    feFunctionArguments args;
    std::vector<double> errorCoeff(4, 0.);
    std::vector<MetricTensor> MX_interpolation(_nVerticesPerElmOnBackmesh);
    std::vector<MetricTensor> MY_interpolation(_nVerticesPerElmOnBackmesh);

#if defined(HAVE_OMP)
#pragma omp for reduction(+ : Area, E_DROIT, E_CURVE, COMP)
#endif
    for(size_t iElm = 0; iElm < nElm; iElm++) {

      double errorLinearOnElm = 0.;
      double errorCurvedOnElm = 0.;

      if(!withAnalyticMetric) {
        for(int i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
          tags[i] = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + i];
          MX_interpolation[i].assignMatrixFrom(_dmmud_dx.at(tags[i]));
          MY_interpolation[i].assignMatrixFrom(_dmmud_dy.at(tags[i]));
        }

        logM0.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[0]));
        logM1.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[1]));
        logM2.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[2]));
      }

      for(size_t i = 0; i < nQuad; i++) {
        xsi[0] = localCoord[3 * i + 0];
        xsi[1] = localCoord[3 * i + 1];

        // Coordonnees physiques
        args.pos[0] = pts[iElm * 3 * nQuad + 3 * i + 0];
        args.pos[1] = pts[iElm * 3 * nQuad + 3 * i + 1];
        args.pos[2] = 0.;

        _options.firstDerivatives->eval(args, DU_EXACT);
        _options.secondDerivatives->eval(args, D2U_EXACT);
        _options.thirdDerivatives->eval(args, D3U_EXACT);

        G[0] = DU_EXACT[0];
        G[1] = DU_EXACT[1];
        
        H(0,0) = D2U_EXACT[0];
        H(0,1) = D2U_EXACT[1];
        H(1,0) = D2U_EXACT[2];
        H(1,1) = D2U_EXACT[3];

        if(withAnalyticMetric) {

          // Evaluate analytic metric and derivatives of M^(-1/2)
          _options.analyticMetric->eval(args, M);
          _options.analytic_dMmud_dx->eval(args, dMmud_dx);
          _options.analytic_dMmud_dy->eval(args, dMmud_dy);

        } else {

          // Interpolate the metric using log-euclidian interpolation
          if(_nVerticesPerElmOnBackmesh == 3) {
            logEuclidianP1Interpolation(xsi, logM0, logM1, logM2, M);
          }

          if(_nVerticesPerElmOnBackmesh == 6) {
            logM3.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[3]));
            logM4.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[4]));
            logM5.assignMatrixFrom(_logMetricTensorAtNodetags.at(tags[5]));
            logEuclidianP2Interpolation(xsi, logM0, logM1, logM2, logM3, logM4, logM5, M);
          }

          // Interpolate the derivagives of M^(-1/2) using the shape functions
          // Derivatives are not SPD, so cant use log-euclidian interpolation
          if(_nVerticesPerElmOnBackmesh == 3) {
            dMmud_dx.assignMatrixFrom(MX_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 0]
                                    + MX_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 1]
                                    + MX_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 2]);
            dMmud_dy.assignMatrixFrom(MY_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 0]
                                    + MY_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 1]
                                    + MY_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 2]);
          } else if(_nVerticesPerElmOnBackmesh == 6) {
            dMmud_dx.assignMatrixFrom(MX_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 0]
                                    + MX_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 1]
                                    + MX_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 2]
                                    + MX_interpolation[3] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 3]
                                    + MX_interpolation[4] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 4]
                                    + MX_interpolation[5] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 5]);
            dMmud_dy.assignMatrixFrom(MY_interpolation[0] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 0]
                                    + MY_interpolation[1] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 1]
                                    + MY_interpolation[2] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 2]
                                    + MY_interpolation[3] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 3]
                                    + MY_interpolation[4] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 4]
                                    + MY_interpolation[5] * basisFunctions[_nVerticesPerElmOnBackmesh * i + 5]);
          } else {
            feInfo("Number of vertices per element on backmesh not supported. Expected P1 or P2 triangle.");
            exit(-1);
          }
        }

        Mmud.assignMatrixFrom(M.pow(-0.5));

        double eloc = -1., eloc_curved = -1.;

        if(computeError) {
          if(_options.polynomialDegree == 1) {
            absH.assignMatrixFrom(H.absoluteValueEigen());
            Eloc.assignMatrixFrom(Mmud * absH * Mmud);

            double contraction = 0.;
            for(int mm = 0; mm < 2; ++mm) {
              for(int ll = 0; ll < 2; ++ll) {
                for(int kk = 0; kk < 2; ++kk) {
                  if(mm == 0) {
                    contraction += G[kk] * Mmud(mm,ll) * dMmud_dx(ll,kk);
                  } else {
                    contraction += G[kk] * Mmud(mm,ll) * dMmud_dy(ll,kk);
                  }
                }
              }
            }

            // Continuous P1 error for straight triangles: 1/8 * trace(Mmud * |H| * Mmud)
            eloc = 1./8. * (Eloc(0,0) + Eloc(1,1));

            // Continuous P1 error for P2 triangles: c * |trace(Mmud * H * Mmud) + (dM^(-1/2) : M^(-1/2)) . gradf|
            Eloc_curved.assignMatrixFrom(Mmud * H * Mmud);
            eloc_curved = 1./8. * fabs(Eloc_curved(0,0) + Eloc_curved(1,1) + contraction);
          } else if(_options.polynomialDegree == 2) {

            // Error term in the reference space for quadratic interpolation on quadratic element
            double resLinear[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
            double resCurved[2][2][2] = {{{0., 0.}, {0., 0.}}, {{0., 0.}, {0., 0.}}};
            double Harray[2][2] = {{D2U_EXACT[0], D2U_EXACT[1]},
                              {D2U_EXACT[2], D2U_EXACT[3]}};
            double C[2][2][2] = {{{D3U_EXACT[0], D3U_EXACT[1]}, {D3U_EXACT[2], D3U_EXACT[3]}},
                                 {{D3U_EXACT[4], D3U_EXACT[5]}, {D3U_EXACT[6], D3U_EXACT[7]}}};

            // For quadratic interpolation: Compute D^3 (f o F(M))_ijk = C_mno * M^(-1/2)_im * M^(-1/2)_jn * M^(-1/2)_ko + ..... LOWER DERIVATIVES
            for(int I = 0; I < 2; ++I) {
              for(int J = 0; J < 2; ++J) {
                for(int K = 0; K < 2; ++K) {
                  for(int m = 0; m < 2; ++m) {
                    for(int n = 0; n < 2; ++n) {
                      for(int o = 0; o < 2; ++o) {

                        // The term for linear geometry
                        resLinear[I][J][K] += C[m][n][o] * Mmud(I,m) * Mmud(J,n) * Mmud(K,o);
                        resCurved[I][J][K] += C[m][n][o] * Mmud(I,m) * Mmud(J,n) * Mmud(K,o);

                        if(_options.curvedMetricReferenceSpaceOptions.enableCurvature) {

                          if(o == 0) {
                            resCurved[I][J][K] += (Harray[m][n] * dMmud_dx(m,I) * Mmud(o,J) * Mmud(n,K)
                                                 + Harray[m][n] * dMmud_dx(m,K) * Mmud(o,I) * Mmud(n,J)
                                                 + Harray[m][n] * dMmud_dx(m,J) * Mmud(o,K) * Mmud(n,I));
                          } else {
                            resCurved[I][J][K] += (Harray[m][n] * dMmud_dy(m,I) * Mmud(o,J) * Mmud(n,K)
                                                 + Harray[m][n] * dMmud_dy(m,K) * Mmud(o,I) * Mmud(n,J)
                                                 + Harray[m][n] * dMmud_dy(m,J) * Mmud(o,K) * Mmud(n,I));
                          }
                        }

                      }
                    }
                  }
                }
              }
            }

            // Continuous error for linear triangles
            errorCoeff[0] =  resLinear[0][0][0];
            errorCoeff[1] =  resLinear[0][0][1] + resLinear[0][1][0] + resLinear[1][0][0];
            errorCoeff[2] =  resLinear[0][1][1] + resLinear[1][0][1] + resLinear[1][1][0];
            errorCoeff[3] =  resLinear[1][1][1];

            bool success = computeAnalyticMetricP2ForLpNorm(errorCoeff, Qtri, 100.);
            if(!success) {
              feInfo("Could not compute upper bound for quadratic error term in reference space!");
              exit(-1);
            }

            eloc = 1./8. * fabs(Qtri(0,0) + Qtri(1,1));

            // Continuous error for quadratic triangles
            errorCoeff[0] =  resCurved[0][0][0];
            errorCoeff[1] =  resCurved[0][0][1] + resCurved[0][1][0] + resCurved[1][0][0];
            errorCoeff[2] =  resCurved[0][1][1] + resCurved[1][0][1] + resCurved[1][1][0];
            errorCoeff[3] =  resCurved[1][1][1];

            success = computeAnalyticMetricP2ForLpNorm(errorCoeff, Qtri, 100.);
            if(!success) {
              feInfo("Could not compute upper bound for quadratic error term in reference space!");
              exit(-1);
            }

            eloc_curved = 1./8. * fabs(Qtri(0,0) + Qtri(1,1));

          } else {
            feInfo("Implement continuous error model pour polynomials of degree %d", _options.polynomialDegree);
            exit(-1);
          }
        }

        Area += weights[i] * detJac[iElm * nQuad + i];
        E_DROIT += weights[i] * detJac[iElm * nQuad + i] * eloc;
        E_CURVE += weights[i] * detJac[iElm * nQuad + i] * eloc_curved;
        COMP += weights[i] * detJac[iElm * nQuad + i] * sqrt(M.determinant());

        errorLinearOnElm += weights[i] * detJac[iElm * nQuad + i] * eloc;
        errorCurvedOnElm += weights[i] * detJac[iElm * nQuad + i] * eloc_curved;
      }

      #pragma omp critical
      {
        errorLinearOnElements[iElm] = errorLinearOnElm;
        errorCurvedOnElements[iElm] = errorCurvedOnElm;
      }

      // Plot error on element
      if(plotError && _nVerticesPerElmOnBackmesh == 3) {
        fprintf(fileLinear, "ST(%+-1.6e,%+-1.6e,0.,%+-1.6e,%+-1.6e,0.,%+-1.6e,%+-1.6e,0.){%+-1.6e,%+-1.6e,%+-1.6e};\n",
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 0] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 0] - 1) + 1],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 1] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 1] - 1) + 1],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 2] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 2] - 1) + 1],
          errorLinearOnElm, errorLinearOnElm, errorLinearOnElm);
        fprintf(fileCurved, "ST(%+-1.6e,%+-1.6e,0.,%+-1.6e,%+-1.6e,0.,%+-1.6e,%+-1.6e,0.){%+-1.6e,%+-1.6e,%+-1.6e};\n",
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 0] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 0] - 1) + 1],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 1] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 1] - 1) + 1],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 2] - 1) + 0],
          coord[3 * (elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + 2] - 1) + 1],
          errorCurvedOnElm, errorCurvedOnElm, errorCurvedOnElm);
      }

    } // for elm
  } // omp parallel

  // E_DROIT = sqrt(E_DROIT);
  // E_CURVE = sqrt(E_CURVE);

  // Ignore the elements with the top 5% error
  std::sort(errorLinearOnElements.begin(), errorLinearOnElements.end());
  size_t numElmToIgnore = ceil(errorLinearOnElements.size() * 0.05);
  double errorLinearTrimmed = 0.;
  for(size_t i = 0; i < errorLinearOnElements.size() - numElmToIgnore; ++i)
    errorLinearTrimmed += errorLinearOnElements[i];
  // errorLinearTrimmed = sqrt(errorLinearTrimmed);
  std::sort(errorCurvedOnElements.begin(), errorCurvedOnElements.end());
  numElmToIgnore = ceil(errorCurvedOnElements.size() * 0.05);
  double errorCurvedTrimmed = 0.;
  for(size_t i = 0; i < errorCurvedOnElements.size() - numElmToIgnore; ++i)
    errorCurvedTrimmed += errorCurvedOnElements[i];
  // errorCurvedTrimmed = sqrt(errorCurvedTrimmed);

  feInfoCond(FE_VERBOSE > 0, "Integral = %+-1.16e", Area);
  feInfoCond(FE_VERBOSE > 0, "E_DROIT = %+-1.16e - Without top 5pc = %+-1.16e", E_DROIT, errorLinearTrimmed);
  feInfoCond(FE_VERBOSE > 0, "E_CURVE = %+-1.16e - Without top 5pc = %+-1.16e --- E_DROIT/E_CURVE = %f - Without top 5pc = %f", 
    E_CURVE, errorCurvedTrimmed, E_DROIT/E_CURVE, errorLinearTrimmed/errorCurvedTrimmed);
  feInfoCond(FE_VERBOSE > 0, "COMP = %+-1.16e", COMP);

  // errorStraight = E_DROIT;
  // errorCurved = E_CURVE;
  errorStraight = errorLinearTrimmed;
  errorCurved = errorCurvedTrimmed;
  complexity = COMP;

  if(plotError) {
    fprintf(fileLinear, "};"); fclose(fileLinear);
    fprintf(fileCurved, "};"); fclose(fileCurved);
  }

#else
  UNUSED(withAnalyticMetric, computeError, plotError, errorFileName, errorStraight, errorCurved, complexity);
#endif
}

// Compute the metric tensor field from recovered derivatives,
// scale it to obtain roughly N vertices and write it as a Gmsh view
// to adapt with MMG
feStatus feMetric::computeMetrics()
{
#if defined(HAVE_GMSH)

  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Computing metrics on Gmsh model %s",
             _options.modelForMetric.c_str());
  gmsh::model::setCurrent(_options.modelForMetric);

  std::vector<std::size_t> nodeTags;
  std::vector<double> coord;
  feStatus s = createVertex2NodeMap(nodeTags, coord);
  if(s != FE_STATUS_OK) return s;

  double N = (double)_options.nTargetVertices;
  double p = _options.LpNorm;
  double n = 2.; // Space dimension
  double deg = (double)_options.polynomialDegree;

  bool scaleMetricField = true;

  //////////////////////////////////////////////////////////////
  // To test convergence of the metric interpolation (for the thesis?):
  // interpolationTest(nodeTags, coord);
  // return FE_STATUS_OK;
  //////////////////////////////////////////////////////////////

  // Compute the raw metric field (no scaling)
  switch(_options.method) {

    case adaptationMethod::ANALYTIC_METRIC:
      s = setAnalyticMetric(nodeTags, coord);
      scaleMetricField = false;
      break;

    case adaptationMethod::ISO_P1: {
      bool isotropic = true;
      s = computeMetricsP1(nodeTags, coord, isotropic);
      break;
    }

    case adaptationMethod::ANISO_P1:
      s = computeMetricsP1(nodeTags, coord);
      break;

    case adaptationMethod::ISO_PN: {
      bool isotropic = true;
      s = computeMetricsPn(nodeTags, coord, isotropic);
      break;
    }

    case adaptationMethod::ANISO_P2:
      s = computeMetricsP2(nodeTags, coord);
      break;

    case adaptationMethod::ANISO_P1_REFERENCE_SPACE:
      {
        int maxIter = 10;
        for(int iter = 0; iter < maxIter; ++iter) {
          // Compute Q_triangle in the reference space, then Q = M^(1/2) * Q_triangle * M^(1/2)
          bool setToIdentity = (iter == 0);
          s = computeMetricsP1_referenceSpace(nodeTags, coord, setToIdentity);
          if(s != FE_STATUS_OK) return s;
          std::string name1 = "Q_iter" + std::to_string(iter) + ".pos";
          std::string name2 = "M_iter" + std::to_string(iter) + ".pos";
          drawEllipsoids(name1, _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
          // Get M from Q
          // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
          // drawEllipsoids(name2, _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        }
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_final.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);

        // // Compute in the physical space to compare
        // s = computeMetricsP2(nodeTags, coord);
        // drawEllipsoids("Q_physical.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_physical.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        exit(-1);
        break;
      }

    case adaptationMethod::ANISO_P2_REFERENCE_SPACE:
      {
        int maxIter = 10;
        for(int iter = 0; iter < maxIter; ++iter) {
          // Compute Q_triangle in the reference space, then Q = M^(1/2) * Q_triangle * M^(1/2)
          bool setToIdentity = (iter == 0);
          s = computeMetricsP2_referenceSpace(nodeTags, coord, setToIdentity);
          if(s != FE_STATUS_OK) return s;
          std::string name1 = "Q_iter" + std::to_string(iter) + ".pos";
          std::string name2 = "M_iter" + std::to_string(iter) + ".pos";
          // drawEllipsoids(name1, _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
          // Get M from Q
          // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
          // drawEllipsoids(name2, _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        }
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_final.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);

        // Compute in the physical space to compare
        s = computeMetricsP2(nodeTags, coord);
        // drawEllipsoids("Q_physical.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_physical.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
        exit(-1);
        break;
      }

    case adaptationMethod::ANISO_PN:
        s = computeMetricsPn(nodeTags, coord);
      break;

    case adaptationMethod::GOAL_ORIENTED_ANISO_P1:
      s = computeMetricsGoalOrientedP1(nodeTags, coord);
      break;

    case adaptationMethod::CURVED_LS:
      s = computeMetricsCurvedLogSimplex(nodeTags, coord);
      break;

    case adaptationMethod::CURVED_LS_INDUCED_DIRECTIONS: {
      bool useInducedDirections = true;
      s = computeMetricsCurvedLogSimplex(nodeTags, coord, useInducedDirections);
      break;
    }

    case adaptationMethod::CURVED_ISOLINES:
      s = computeMetricsCurvedIsolines(nodeTags, coord);
      scaleMetricField = _options.scaleIsolineMetric;
      break;

    case adaptationMethod::CURVED_REFERENCE_SPACE:
      s = computeMetricsCurvedReferenceSpace(nodeTags, coord);
      break;

    case adaptationMethod::CURVED_GRAPH_SURFACE:
      s = computeMetricsP2_forGraphSurface(nodeTags, coord);
      scaleMetricField = false;
      break;
  }

  if(s != FE_STATUS_OK) return s;

  // Reconstruct the underlying manifold before scaling
  // reconstructManifoldFromMetric(nodeTags, coord, _metricTensorAtNodetags, _sequentialTag2nodeTagVec, _recoveredFields[0]);

#if !defined(ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE)
  // Apply scaling and gradation only when not testing for metric interpolation convergence
  if(_options.debug) drawEllipsoids("rawMetrics.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);

  if(scaleMetricField)
  {
    double exponentInIntegral, exponentForDeterminant;

    if(_options.targetNorm == Norm::Lp) {
      exponentInIntegral = p * (deg + 1.0) / (2.* (p * (deg + 1.0) + n));
      exponentForDeterminant = -1. / (p * (deg + 1.) + n);
      feInfo("Exponents for scaling for L%2.1f norm = %+-1.4e - %+-1.4e", p, exponentForDeterminant, exponentInIntegral);
    } else {
      exponentInIntegral = p * deg / (2.* (p * deg + n));
      exponentForDeterminant = -1. / (p * deg + n);
      feInfo("Exponents for scaling for H1 seminorm = %+-1.4e - %+-1.4e", exponentForDeterminant, exponentInIntegral);
    }

    // Scale the metric field
    feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Targetting %f vertices in L%f norm", N, p);
    metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);

    // Bound the eigenvalues according to hMin and hMax
    for(auto &pair : _metricTensorAtNodetags) {
      pair.second = pair.second.boundEigenvalues(_lambdaMin, _lambdaMax);
    }

    if(_options.debug) {
      drawEllipsoids("metricsAfterScaling.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
      if(_options.method == adaptationMethod::ANISO_P2)
        drawEllipsoids("metricsAfterScaling_p2_analytic.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
      if(_options.method == adaptationMethod::ANISO_PN)
        drawEllipsoids("metricsAfterScaling_logSimplex.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
    }
  }

  // Apply gradation
  if(_options.enableGradation) {
    newGradation(nodeTags, coord, _options.gradation, _metricTensorAtNodetags);
  }

  if(_options.debug) {
    drawEllipsoids("metricsAfterGradation.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
  }
#endif

  // Precompute all metric logarithms
  _logMetricTensorAtNodetags_eigen.resize(_metricTensorAtNodetags.size());
  for(auto &pair : _metricTensorAtNodetags) {
    _logMetricTensorAtNodetags[pair.first] = pair.second.log();
    _logMetricTensorAtNodetags_eigen[pair.first - 1] = convert2eigenMatrix(pair.second.log());
  }

  // Write the metric field as a view in the Gmsh model
  writeMetricField(nodeTags, coord);

#endif
  return FE_STATUS_OK;
}