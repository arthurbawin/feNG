
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
  : _recovery(recovery), _options(metricOptions)
{
  _options.polynomialDegree = _recovery->getDegreeSolution();
}

feMetric::feMetric(std::vector<feNewRecovery *> &recoveredFields, feMetricOptions metricOptions)
  : _recoveredFields(recoveredFields), _options(metricOptions)
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
#endif
  return FE_STATUS_OK;
}

// Scale the metric field to fit N vertices in the final mesh. Interpolation is
// performed on the Gmsh subsitute (gmshModel). Using gmsh api because it's easier.
// FIXME: parallelize
template <class MetricType>
void feMetric::metricScalingFromGmshSubstitute(std::map<int, MetricType> &metrics,
                                               const std::vector<size_t> &nodeTags,
                                               const std::vector<double> &coord,
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

      for(size_t i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
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

    for(size_t i = 0; i < nQuad; i++) {
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

      double uexact = matfun(xphys, yphys);

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
    system(cmd.data());
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

#endif
}

// Compute unscaled optimal metric field minimizing interpolation error
// in Lp norm for P1 elements according to Alauzet & Loseille
feStatus feMetric::computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                                    bool isotropic)
{
  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double x[2], fxx, fxy, fyx, fyy;
    MetricTensor H, Q;
    std::vector<double> d2exact(4, 0.), pos(3, 0.);

    // Compute bounded absolute value of hessian at vertices (at nodetags)
    #if defined(HAVE_OMP)
    #pragma omp for
    #endif
    for(size_t i = 0; i < nodeTags.size(); i++) {
      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      // Evaluate 2nd order derivatives
      fxx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 2, x);
      fxy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 3, x);
      fyx = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 4, x);
      fyy = _recoveredFields[0]->evaluateRecovery(PPR::DERIVATIVE, 5, x);

     if(_options.useAnalyticDerivatives) {
        // Evaluate the analytic high-order derivatives instead
        pos[0] = x[0];
        pos[1] = x[1];
        _options.secondDerivatives->eval(0, pos, d2exact);

        fxx = d2exact[0];
        fxy = d2exact[1];
        fyx = d2exact[2];
        fyy = d2exact[3];
      }

      H(0, 0) = fxx;
      H(0, 1) = (fxy + fyx) / 2.;
      H(1, 0) = (fxy + fyx) / 2.;
      H(1, 1) = fyy;

      if(_options.targetNorm == Norm::Lp) {
        if(isotropic)
          Q = H.boundEigenvaluesOfAbsIsotropic(_lambdaMin, _lambdaMax);
        else {
          // Q = H.absoluteValueEigen();
          Q = H.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
        }
      } else {
        // Target norm to minimize is H1. Use Mirebeau's analytical solution.
        double a = fxx / 2.;
        double b = ((fxy + fyx)/2.) / 2.;
        double c = fyy / 2.;
        if(fabs(a*c - b*b) < 1e-10) {
          feInfo("a*c - b*b = %+-1.4e", a*c - b*b);
          exit(-1);
        }
        double alpha = 100.;
        Q(0, 0) = 4. * (a*a + b*b);
        Q(0, 1) = 4. * (a*b + b*c);
        Q(1, 0) = 4. * (a*b + b*c);
        Q(1, 1) = 4. * (b*b + c*c);
        // Limit the anisotropy (Mirebeau's thesis, Section 3.4.3 equation 3.74
        Q = Q.limitAnisotropy(alpha);
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

  double x[2], fxx, fxy, fyx, fyy, p, px, py;
  MetricTensor H;

  // Test pour la perturbation singuliere
  double epsilon = 1e-2;
  double velocity[2] = {1. - 2. * epsilon, 1. - 2. * epsilon};
  double a = 2. * (1 - epsilon);

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
    px = _recoveredFields[1]->evaluateRecovery(PPR::DERIVATIVE, 0, x);
    py = _recoveredFields[1]->evaluateRecovery(PPR::DERIVATIVE, 1, x);

    double dotProd = px * velocity[0] + py * velocity[1];

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

feStatus feMetric::computeMetricsP2(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
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
    std::vector<double> D3U_EXACT(8, 0.), pos(3, 0.);

// Compute bounded absolute value of upper bound Q at vertices (at nodetags)
#if defined(HAVE_OMP)
#pragma omp for
#endif
    for(size_t i = 0; i < numVertices; i++) {

      x[0] = coord[3 * i + 0];
      x[1] = coord[3 * i + 1];

      // Get the coefficients of the homogeneous error polynomial at vertex
      std::vector<double> &errorCoeff = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];

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
        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = 0.;
        _options.thirdDerivatives->eval(0, pos, D3U_EXACT);
        double uxxx = D3U_EXACT[0];
        double uxxy = D3U_EXACT[1];
        double uxyy = D3U_EXACT[3];
        double uyyy = D3U_EXACT[7];
        errorCoeff[0] = uxxx;
        errorCoeff[1] = 3. * uxxy;
        errorCoeff[2] = 3. * uxyy;
        errorCoeff[3] = uyyy;
      }

      if(_options.targetNorm == Norm::Lp) {
        computeAnalyticMetricP2ForLpNorm(errorCoeff, Q, _options.maxAnalyticEllipseDiameter);

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
        computeAnalyticMetricP2ForH1semiNorm(errorCoeff, Q);
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

static MetricTensor metricAparicioEstrems(double x, double y)
{
  double pi = M_PI;
  double A = 10.;
  double c = sqrt(A * A + 4. * M_PI * M_PI);
  double HMIN = 0.1;
  double b = 10.;

  // Size from the paper
  // double h = HMIN + b * fabs(y);
  // double dhdx = 0.;
  // double dhdy = b * sgn(y);

  double h = HMIN + b * fabs((A * y - cos(2. * M_PI * x)) / c);
  double dhdx = b * (2. * pi * sin(2. * pi * x) * sgn(A * y - cos(2. * pi * x))) / c;
  double dhdy = b * (A * sgn(A * y - cos(2. * pi * x))) / c;

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

// Compute scaled optimal metric field for Pn elements
// according to Coulaud & Loseille using the log-simplex method
feStatus feMetric::computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord,
                                    bool isotropic)
{
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
             "Computing metric tensors (ANISO_PN) for %d vertices...", nodeTags.size());

  if(_options.polynomialDegree > 3) {
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
  int cnter = 0;

#if defined(HAVE_OMP)
#pragma omp parallel shared(cnter, OK)
#endif
  {
    double x[2];
    MetricTensor Q;
    std::vector<double> D3U_EXACT(8, 0.), D4U_EXACT(16, 0.), pos(3, 0.);

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

      // Get the coefficients of the homogeneous error polynomial at vertex
      std::vector<double> &errorCoeff = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];

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
      }

      if(_options.useAnalyticDerivatives) {
        pos[0] = x[0];
        pos[1] = x[1];
        pos[2] = 0.;
        // Evaluate the analytic high-order derivatives instead
        if(_options.polynomialDegree == 2) {
          
          _options.thirdDerivatives->eval(0, pos, D3U_EXACT);

          double uxxx = D3U_EXACT[0];
          double uxxy = D3U_EXACT[1];
          double uxyy = D3U_EXACT[3];
          double uyyy = D3U_EXACT[7];

          errorCoeff[0] = uxxx;
          errorCoeff[1] = uxxy;
          errorCoeff[2] = uxyy;
          errorCoeff[3] = uyyy;

        } else if(_options.polynomialDegree == 3) {

          _options.fourthDerivatives->eval(0, pos, D4U_EXACT);

          errorCoeff[0] = D4U_EXACT[0];
          errorCoeff[1] = D4U_EXACT[1];
          errorCoeff[2] = D4U_EXACT[3];
          errorCoeff[3] = D4U_EXACT[7];
          errorCoeff[4] = D4U_EXACT[15];

        } else {
          feInfo("TEMPORARY ERROR MESSAGE: CANNOT COMPUTE ANALYTIC DERIVATIVES FOR THIS POLYNOMIAL DEGREE");
          exit(-1);
        }
      }

      bool res = false;

      if(_options.targetNorm == Norm::Lp) {
        // Optimal upper bound for |pi(z)| <= 1
        res = computeMetricLogSimplexStraight(x, errorCoeff, _options.polynomialDegree, nTheta,
                                                 maxIter, tol, Q, numIter, myLP);
      } else {
        // Optimal upper bound for |grad(pi)(z)| <= 1
        // if(_options.polynomialDegree == 3) {

        // } else {
          OK = false;
          feErrorMsg(FE_STATUS_ERROR, "Implement optimal metric for H1 seminorm for P4+");
        // }
        res = computeMetricLogSimplexStraight(x, errorCoeff, _options.polynomialDegree, nTheta,
                                                 maxIter, tol, Q, numIter, myLP);
      }

      // If metric computation failed, try again by doubling the number of constraints points in
      // each quadrant (x8 total)
      int retry = 0, maxTry = 5;
      while(!res && retry < maxTry) {
        nTheta *= 2;
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
                   "Could not compute metric. Trying again with %d constraints.", nTheta);
        setUpLinearProblem(myLP, _options, nTheta, true);
        res = computeMetricLogSimplexStraight(x, errorCoeff, _options.polynomialDegree, nTheta,
                                              maxIter, tol, Q, numIter, myLP);
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
        OK = false;
        feErrorMsg(FE_STATUS_ERROR,
                   "Could not compute a metric at (%+-1.5e - %+-1.5e) (vertex %d/%d)", x[0], x[1],
                   i, nodeTags.size());
      }
    }
  }

  if(!OK) {
    return feErrorMsg(FE_STATUS_ERROR, "Could not compute at least one metric tensor");
  }

  return FE_STATUS_OK;
#else
  return feErrorMsg(
    FE_STATUS_ERROR,
    "SoPlex (optimization library) is required to compute metric tensors for Pn adaptation");
#endif
}

thread_local std::vector<double> POS(3, 0.), DU(2, 0.), D2U(4, 0.), D3U(8, 0.);

void feMetric::computeDirectionFieldFromGradient(const double x[2], const int vertex, double grad[2], const double tol)
{
  double gx, gy;

  if(_options.useAnalyticDerivatives) {
    // Evaluate the analytic gradient
    POS[0] = x[0];
    POS[1] = x[1];
    POS[2] = 0.;
    _options.firstDerivatives->eval(0, POS, DU);
    gx = DU[0];
    gy = DU[1];
  } else {
    gx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
    gy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);
  }

  double normGrad = sqrt(gx * gx + gy * gy);

  if(normGrad > tol) {
    double theta1 = atan2(gy, gx);
    grad[0] = cos(theta1);
    grad[1] = sin(theta1);
  } else {
    // Gradient norm is too small : directions will be smoothed during gradation
    grad[0] = 1.;
    grad[1] = 0.;
  }
}

inline double maxFabs(const double x, const double threshold)
{
  if(x < 0 && fabs(x) < threshold) return -threshold;
  if(x > 0 && fabs(x) < threshold) return threshold;
  return x;
}

double feMetric::secondDerivativesAlongCurve(const double x[2], const int vertex, const double directionV1[2], const int direction)
{
  double fx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
  double fy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

  double fxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
  double fxy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
                _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) /
               2.;
  double fyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

  double TOL = 1e-12;

  if(_options.useAnalyticDerivatives) {
    // Evaluate the analytic high-order derivatives instead
    POS[0] = x[0];
    POS[1] = x[1];
    POS[2] = 0.;
    _options.firstDerivatives->eval(0, POS, DU);
    fx = maxFabs(DU[0], TOL);
    fy = maxFabs(DU[1], TOL);
    _options.secondDerivatives->eval(0, POS, D2U);
    fxx = maxFabs(D2U[0], TOL);
    fxy = maxFabs(D2U[1], TOL);
    fyy = maxFabs(D2U[3], TOL);
  }

  double kappa1 = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappa2 = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappa[2] = {kappa1, kappa2};

  double C = directionV1[0];
  double S = directionV1[1];
  double g[2], gOrth[2];

  if(direction == 0) {
    g[0] = C;
    g[1] = S;
    gOrth[0] = -S;
    gOrth[1] = C;
  } else {
    gOrth[0] = C;
    gOrth[1] = S;
    g[0] = -S;
    g[1] = C;
  }

  double d2 = g[0] * g[0] * fxx + 2. * g[0] * g[1] * fxy + g[1] * g[1] * fyy;
  // double d1 = 2. * kappa[direction] * (gOrth[0] * fx + gOrth[1] * fy);
  double d1 = kappa[direction] * (g[0] * fx + g[1] * fy);
  feInfo("dir = %d - k = %f - 1/r = %f - %f - %f - %1.4e", direction, kappa[direction], 1./sqrt(x[0]*x[0] + x[1]*x[1]), d2, d1, d2+d1);

  return d2 + d1;
}

thread_local std::set<double> ROOTS_L_LINEAR;
thread_local Eigen::VectorXd COEFF_L_LINEAR = Eigen::VectorXd::Zero(5);

double feMetric::solveSizePolynomialLinear(const double x[2], const int vertex, const double directionV1[2], const int direction, const double targetError)
{
  double fx, fy, fxx, fxy, fyy, TOL = 1e-12;

  if(_options.useAnalyticDerivatives) {
    // Evaluate the analytic high-order derivatives instead
    POS[0] = x[0];
    POS[1] = x[1];
    POS[2] = 0.;
    _options.firstDerivatives->eval(0, POS, DU);
    fx = maxFabs(DU[0], TOL);
    fy = maxFabs(DU[1], TOL);
    _options.secondDerivatives->eval(0, POS, D2U);
    fxx = maxFabs(D2U[0], TOL);
    fxy = maxFabs(D2U[1], TOL);
    fyy = maxFabs(D2U[3], TOL);
  } else {
    fx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
    fy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

    fxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
    fxy = 0.5 * (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex));
    fyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);
  }

  // Compute directions and curvatures
  double ax, ay, bx, by, kappa;

  // (C,S) is the direction of the gradient at x
  // double C = directionV1[0];
  // double S = directionV1[1];

  double normGrad = fmax(TOL, sqrt(fx*fx+fy*fy));

  double  xIso[2] = {-fy,  fy};
  double  yIso[2] = { fx, -fx};
  double xGrad[2] = { fx, -fx};
  double yGrad[2] = { fy, -fy};
  double  kappaIso = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappaGrad = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappaVec[2] = {kappaIso, kappaGrad};
  kappa = kappaVec[direction];

  // // 0 is iso, 1 is gradient
  // if(direction == 0) {
  //   kappa = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  //   ax = -fy/normGrad;
  //   ay =  fx/normGrad;
  //   bx =  fx/normGrad * 0.5 * kappa;
  //   by =  fy/normGrad * 0.5 * kappa;
  // } else {
  //   kappa = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
  //   ax =  fx/normGrad;
  //   ay =  fy/normGrad;
  //   bx = -fy/normGrad * 0.5 * kappa;
  //   by =  fx/normGrad * 0.5 * kappa;
  // }

  double L = 1e10;

  // Forward and backward along iso or grad
  for(int i = 0; i < 1; ++i) {

    if(direction == 0) {
      ax =  xIso[i]/normGrad;
      ay =  yIso[i]/normGrad;
      bx = xGrad[i]/normGrad * 0.5 * kappa;
      by = yGrad[i]/normGrad * 0.5 * kappa;
    } else {
      ax = xGrad[i]/normGrad;
      ay = yGrad[i]/normGrad;
      bx =  xIso[i]/normGrad * 0.5 * kappa;
      by =  yIso[i]/normGrad * 0.5 * kappa;
    }

    double coeffL4 = (fxx * bx * bx + 2. * fxy * bx*by + fyy*by*by) / 3.;
    double coeffL3 = (ax*bx*fxx + fxy * (ax*by + ay*bx) + ay*by*fyy) * 2./3.;
    double coeffL2 = (0.5 * (fxx*ax*ax + fyy*ay*ay) + fxy*ax*ay + bx*fx + by*fy);

    // Solve p(L) = error
    COEFF_L_LINEAR(0) = coeffL4;
    COEFF_L_LINEAR(1) = coeffL3;
    COEFF_L_LINEAR(2) = coeffL2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    // feInfo("Polynomial is %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", 
    //   COEFF_L_LINEAR(0),
    //   COEFF_L_LINEAR(1),
    //   COEFF_L_LINEAR(2),
    //   COEFF_L_LINEAR(3),
    //   COEFF_L_LINEAR(4));
    ROOTS_L_LINEAR = RootFinder::solvePolynomial(COEFF_L_LINEAR, 0., INFINITY, 1e-4);
    for(auto it = ROOTS_L_LINEAR.begin(); it != ROOTS_L_LINEAR.end(); it++) {
      L = fmin(L, fabs(*it));
    }

    // Solve -p(L) = error
    COEFF_L_LINEAR(0) = -coeffL4;
    COEFF_L_LINEAR(1) = -coeffL3;
    COEFF_L_LINEAR(2) = -coeffL2;
    COEFF_L_LINEAR(3) = 0.;
    COEFF_L_LINEAR(4) = -targetError;
    // feInfo("Polynomial is %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", 
    //   COEFF_L_LINEAR(0),
    //   COEFF_L_LINEAR(1),
    //   COEFF_L_LINEAR(2),
    //   COEFF_L_LINEAR(3),
    //   COEFF_L_LINEAR(4));
    ROOTS_L_LINEAR = RootFinder::solvePolynomial(COEFF_L_LINEAR, 0., INFINITY, 1e-4);
    for(auto it = ROOTS_L_LINEAR.begin(); it != ROOTS_L_LINEAR.end(); it++) {
      L = fmin(L, fabs(*it));
    }
  }

  return L;
}

thread_local std::set<double> ROOTS_L_QUADRATIC;
thread_local Eigen::VectorXd COEFF_L_QUADRATIC = Eigen::VectorXd::Zero(7);

double feMetric::solveSizePolynomialQuadratic(const double x[2], const int vertex, const double directionV1[2], const int direction, const double targetError)
{
  double fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy, TOL = 1e-12;

  if(_options.useAnalyticDerivatives) {
    // Evaluate the analytic high-order derivatives instead
    POS[0] = x[0];
    POS[1] = x[1];
    POS[2] = 0.;
    _options.firstDerivatives->eval(0, POS, DU);
    fx = maxFabs(DU[0], TOL);
    fy = maxFabs(DU[1], TOL);
    _options.secondDerivatives->eval(0, POS, D2U);
    fxx = maxFabs(D2U[0], TOL);
    fxy = maxFabs(D2U[1], TOL);
    fyy = maxFabs(D2U[3], TOL);
    _options.thirdDerivatives->eval(0, POS, D3U);
    fxxx = maxFabs(D3U[0], TOL);
    fxxy = maxFabs(D3U[1], TOL);
    fxyy = maxFabs(D3U[3], TOL);
    fyyy = maxFabs(D3U[7], TOL);
  } else {
    fx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
    fy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

    fxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
    fxy = 0.5 * (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex));
    fyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

    fxxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
    fxxy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex) +
            _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex) +
            _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex)) / 3.;
    fxyy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex) +
            _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex) +
            _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex)) / 3.;
    fyyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);
  }

  // Compute directions and curvatures
  double ax, ay, bx, by, kappa;

  // (C,S) is the direction of the gradient at x
  // double C = directionV1[0];
  // double S = directionV1[1];

  double normGrad = fmax(TOL, sqrt(fx*fx+fy*fy));

  double  xIso[2] = {-fy,  fy};
  double  yIso[2] = { fx, -fx};
  double xGrad[2] = { fx, -fx};
  double yGrad[2] = { fy, -fy};
  double  kappaIso = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappaGrad = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappaVec[2] = {kappaIso, kappaGrad};
  kappa = kappaVec[direction];

  double L = 1e10;

  // Forward and backward along iso or grad
  for(int i = 0; i < 1; ++i) {

    if(direction == 0) {
      ax =  xIso[i]/normGrad;
      ay =  yIso[i]/normGrad;
      bx = xGrad[i]/normGrad * 0.5 * kappa;
      by = yGrad[i]/normGrad * 0.5 * kappa;
    } else {
      ax = xGrad[i]/normGrad;
      ay = yGrad[i]/normGrad;
      bx =  xIso[i]/normGrad * 0.5 * kappa;
      by =  yIso[i]/normGrad * 0.5 * kappa;
    }

    double coeffL6 = (fxxx * bx*bx*bx + 3. * fxxy * bx*bx*by + 3. * fxyy * bx*by*by + fyyy*by*by*by) * 2./15.;
    double coeffL5 = (fxxx * bx*bx*ax +      fxxy * bx*bx*ay +      fxyy * ax*by*by + fyyy*ay*by*by) * 2./5.
                    + (fxxy * ax*bx*by + fxyy * ay*bx*by) * 4./5.;
    double coeffL4 = (fxxx * ax*ax*bx /2. + fxxy * ax*ax*by/2. + fxxy*ax*ay*bx + fxyy*ax*ay*by + fxyy*ay*ay*bx/2. + fyyy*ay*ay*by/2. + fxx*bx*bx + 2.*fxy*bx*by + fyy*by*by);
    double coeffL3 = (fxxx*ax*ax*ax/3. + fyyy*ay*ay*ay/3. + 2. * ax * (bx*fxx + by*fxy) + 2.*ay*(bx*fxy + by*fyy) + ax*ax*ay*fxxy + ax*ay*ay*fyy);

    // Solve p(L) = error
    COEFF_L_QUADRATIC(0) = coeffL6;
    COEFF_L_QUADRATIC(1) = coeffL5;
    COEFF_L_QUADRATIC(2) = coeffL4;
    COEFF_L_QUADRATIC(3) = coeffL3;
    COEFF_L_QUADRATIC(4) = 0.;
    COEFF_L_QUADRATIC(5) = 0.;
    COEFF_L_QUADRATIC(6) = -targetError;
    // feInfo("Polynomial is %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", 
    //   COEFF_L_QUADRATIC(0),
    //   COEFF_L_QUADRATIC(1),
    //   COEFF_L_QUADRATIC(2),
    //   COEFF_L_QUADRATIC(3),
    //   COEFF_L_QUADRATIC(4));
    ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, 0., INFINITY, 1e-4);
    for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
      L = fmin(L, fabs(*it));
    }

    // Solve -p(L) = error
    COEFF_L_QUADRATIC(0) = -coeffL6;
    COEFF_L_QUADRATIC(1) = -coeffL5;
    COEFF_L_QUADRATIC(2) = -coeffL4;
    COEFF_L_QUADRATIC(3) = -coeffL3;
    COEFF_L_QUADRATIC(4) = 0.;
    COEFF_L_QUADRATIC(5) = 0.;
    COEFF_L_QUADRATIC(6) = -targetError;
    // feInfo("Polynomial is %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e - %+-1.4e", 
    //   COEFF_L_QUADRATIC(0),
    //   COEFF_L_QUADRATIC(1),
    //   COEFF_L_QUADRATIC(2),
    //   COEFF_L_QUADRATIC(3),
    //   COEFF_L_QUADRATIC(4));
    ROOTS_L_QUADRATIC = RootFinder::solvePolynomial(COEFF_L_QUADRATIC, 0., INFINITY, 1e-4);
    for(auto it = ROOTS_L_QUADRATIC.begin(); it != ROOTS_L_QUADRATIC.end(); it++) {
      L = fmin(L, fabs(*it));
    }
  }

  return L;
}

double feMetric::thirdDerivativesAlongCurve(const double x[2], const int vertex, const double directionV1[2], const int direction)
{
  double fx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 0, vertex);
  double fy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 1, vertex);

  double fxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 2, vertex);
  double fxy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 3, vertex) +
                _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 4, vertex)) /
               2.;
  double fyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 5, vertex);

  double fxxx = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 6, vertex);
  double fxxy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 7, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 8, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 10, vertex)) /
                3.;
  double fxyy = (_recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 9, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 11, vertex) +
                 _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 12, vertex)) /
                3.;
  double fyyy = _recoveredFields[0]->evaluateRecoveryAtVertex(PPR::DERIVATIVE, 13, vertex);

  double TOL = 1e-12;

  if(_options.useAnalyticDerivatives) {
    // Evaluate the analytic high-order derivatives instead
    POS[0] = x[0];
    POS[1] = x[1];
    POS[2] = 0.;
    _options.firstDerivatives->eval(0, POS, DU);
    fx = maxFabs(DU[0], TOL);
    fy = maxFabs(DU[1], TOL);
    _options.secondDerivatives->eval(0, POS, D2U);
    fxx = maxFabs(D2U[0], TOL);
    fxy = maxFabs(D2U[1], TOL);
    fyy = maxFabs(D2U[3], TOL);
    _options.thirdDerivatives->eval(0, POS, D3U);
    fxxx = maxFabs(D3U[0], TOL);
    fxxy = maxFabs(D3U[1], TOL);
    fxyy = maxFabs(D3U[3], TOL);
    fyyy = maxFabs(D3U[7], TOL);
  }

  double kappa1 = (-fy * fy * fxx + 2.0 * fx * fy * fxy - fx * fx * fyy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappa2 = (fx * fy * (fyy - fxx) + (fx * fx - fy * fy) * fxy) / (pow(fx * fx + fy * fy, 3. / 2.));
  double kappa[2] = {kappa1, kappa2};

  double C = directionV1[0];
  double S = directionV1[1];
  double g[2], gOrth[2];

  if(direction == 0) {
    g[0] = C;
    g[1] = S;
    gOrth[0] = -S;
    gOrth[1] = C;
  } else {
    gOrth[0] = C;
    gOrth[1] = S;
    g[0] = -S;
    g[1] = C;
  }

  double d3 = g[0] * g[0] * g[0] * fxxx + 3. * g[0] * g[0] * g[1] * fxxy +
              3. * g[0] * g[1] * g[1] * fxyy + g[1] * g[1] * g[1] * fyyy;

  double d2 = 3. * kappa[direction] * (g[0] * gOrth[0] * fxx + (g[0] * gOrth[1] + g[1] * gOrth[0]) * fxy + g[1] * gOrth[1] * fyy);

  return d3 + d2;
}

 void feMetric::computeSizeField(const double pos[2], const int vertex, const double directionV1[2], double &hGrad, double &hIso)
{
  switch(_options.polynomialDegree) {
    case 1:
      {
        // // // double derivativeAlongIsoline = fabs(secondDerivativesAlongCurve(pos, vertex, directionV1, 0));
        // // // double derivativeAlongGradient = fabs(secondDerivativesAlongCurve(pos, vertex, directionV1, 1));
        // // // hIso  = sqrt(2.0 * _options.eTargetError / fmax(derivativeAlongIsoline, 1e-10));
        // // // hGrad = sqrt(2.0 * _options.eTargetError / fmax(derivativeAlongGradient, 1e-10));

        // Test for f(x,y) = x^2+y^2 only
        // hGrad = sqrt(_options.eTargetError);
        // // hIso = sqrt(2. * _options.eTargetError * sqrt(pos[0]*pos[0]+pos[1]*pos[1]));
        // hIso = sqrt(sqrt(6.)) * sqrt(_options.eTargetError * sqrt(pos[0]*pos[0]+pos[1]*pos[1]));

        // hIso  = solveSizePolynomialLinear(pos, vertex, directionV1, 0, pow(_options.eTargetError, 2) );
        // hGrad = solveSizePolynomialLinear(pos, vertex, directionV1, 1,     _options.eTargetError     );

        hIso  = solveSizePolynomialLinear(pos, vertex, directionV1, 0, _options.eTargetError);
        hGrad = solveSizePolynomialLinear(pos, vertex, directionV1, 1, _options.eTargetError);

        // if(fabs(hIso_ini - hIso)/fabs(hIso_ini) > 1e-4 || fabs(hGrad_ini - hGrad)/fabs(hGrad_ini) > 1e-4) {
        //   feInfo("At (%+-1.4e,%+-1.4e) : %+-1.4e vs %+-1.4e -- %+-1.4e vs %+-1.4e", pos[0], pos[1], hIso_ini, hIso, hGrad_ini, hGrad);
        //   feErrorMsg(FE_STATUS_ERROR, "Human is dead. Mismatch.");
        //   exit(-1);
        // }

      }
      return;
    case 2:
      {
        // double derivativeAlongGradient = fabs(thirdDerivativesAlongCurve(pos, vertex, directionV1, 0));
        // double derivativeAlongIsoline = fabs(thirdDerivativesAlongCurve(pos, vertex, directionV1, 1));
        // hIso = pow(6.0 * _options.eTargetError / derivativeAlongIsoline, 1./3.);
        // hGrad = pow(6.0 * _options.eTargetError / derivativeAlongGradient, 1./3.);

        hIso  = solveSizePolynomialQuadratic(pos, vertex, directionV1, 0, _options.eTargetError);
        hGrad = solveSizePolynomialQuadratic(pos, vertex, directionV1, 1, _options.eTargetError);

      }
      return;
    default:
      feErrorMsg(FE_STATUS_ERROR, "Cannot compute size field for interpolation of order > 2. Exiting.");
      exit(-1);
  }
}

feStatus feMetric::computeMetricsExtremeSizesOnly(std::vector<std::size_t> &nodeTags,
                                                  std::vector<double> &coord)
{
  size_t numVertices = nodeTags.size();
  // std::vector<double[2]> directionV1(numVertices);
  double x[2], direction[2], C, S;

  FILE *directions, *dirSmoothed;
  directions = fopen("directions.pos", "w");
  dirSmoothed = fopen("directionsSmoothed.pos", "w");
  fprintf(directions, "View \" directions \"{\n");
  fprintf(dirSmoothed, "View \" directionsSmoothed \"{\n");

  std::map<size_t, double> COS;
  std::map<size_t, double> SIN;

  // Compute the direction field
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];
    computeDirectionFieldFromGradient(x, _nodeTag2sequentialTag[nodeTags[i]], direction, 1e-12);
    C = direction[0];
    S = direction[1];
    fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
    fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.);
    COS[nodeTags[i]] = C;
    SIN[nodeTags[i]] = S;
  }

  if(_options.smoothDirectionField) {
    // Smooth the direction field
    int smoothMaxIter = 100;
    double smoothTol = 0.95;
    smoothDirections(COS, SIN,  smoothMaxIter, smoothTol);
  }

  for(size_t i = 0; i < numVertices; i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    C = COS[nodeTags[i]];
    S = SIN[nodeTags[i]];

    if(_options.smoothDirectionField) {
      fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
      fprintf(dirSmoothed, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.);
    } else {
      fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., C, S, 0.);
      fprintf(directions, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.);
    } 

    double hGrad = 1., hIso = 1.;
    direction[0] = C;
    direction[1] = S; 
    computeSizeField(x, _nodeTag2sequentialTag[nodeTags[i]], direction, hGrad, hIso);

    double eigenvalues[2] = {fmax(_lambdaMin, fmin(_lambdaMax, 1. / (hGrad * hGrad))), 
                             fmax(_lambdaMin, fmin(_lambdaMax, 1. / (hIso  * hIso )))};
    double ev1[2] = {C, S};
    double ev2[2] = {-S, C};

    // MetricTensor M(eigenvalues, ev1, ev2);
    // _metricTensorAtNodetags[nodeTags[i]] = M.boundEigenvaluesOfAbs(_lambdaMin, _lambdaMax);
    _metricTensorAtNodetags[nodeTags[i]] = MetricTensor(eigenvalues, ev1, ev2);
  }

  fprintf(directions, "};"); fclose(directions);
  fprintf(dirSmoothed, "};"); fclose(dirSmoothed);

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

      // Get the direction of the gradient
      double grad[2];
      computeDirectionFieldFromGradient(x, _nodeTag2sequentialTag[nodeTags[i]], grad, 1e-8);
      double C = grad[0];
      double S = grad[1];
      // #pragma omp critical
      // {
      //   fprintf(directionFile, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0.,  C, S, 0.);
      //   fprintf(directionFile, "VP(%g,%g,%g){%g,%g,%g};\n", x[0], x[1], 0., -S, C, 0.);
      // }

      bool res =
        computeMetricLogSimplexCurved(_nodeTag2sequentialTag[nodeTags[i]], x, grad,
                                      _recoveredFields[0], Q, maxIter, nTheta, tol, numIter, myLP);

      // If metric computation failed, try again by doubling the number of constraints points in
      // each quadrant (x8 total)
      int retry = 0, maxTry = 5;
      while(!res && retry < maxTry) {
        nTheta *= 2;
        feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE,
                   "Could not compute metric. Trying again with %d constraints.", nTheta);
        setUpLinearProblem(myLP, _options, nTheta, true);
        res = computeMetricLogSimplexCurved(_nodeTag2sequentialTag[nodeTags[i]], x, grad,
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
            double iso[2] = {-grad[1], grad[0]};
            _metricTensorAtNodetags[nodeTags[i]] =
              Q.setEigenvectorsAndBoundEigenvalues(iso, grad, _lambdaMin, _lambdaMax);
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
  return feErrorMsg(
    FE_STATUS_ERROR,
    "SoPlex (optimization library) is required to compute metric tensors for curved P2 adaptation");
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

  //////////////////////////////////////////////////////////////
  // To test convergence of the metric interpolation (for the thesis?):
  // interpolationTest(nodeTags, coord);
  // return FE_STATUS_OK;
  //////////////////////////////////////////////////////////////

  // Compute the raw metric field (no scaling)
  switch(_options.method) {
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

    case adaptationMethod::CURVED_EXTREME_SIZES:
      s = computeMetricsExtremeSizesOnly(nodeTags, coord);
      break;
  }

  if(s != FE_STATUS_OK) return s;

#if !defined(ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE)
  // Apply scaling and gradation only when not testing for metric interpolation convergence
  if(_options.debug) drawEllipsoids("rawMetrics.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);

  if(_options.method != adaptationMethod::CURVED_EXTREME_SIZES) {
    // Scale the metric field
    double N = (double)_options.nTargetVertices;
    double p = _options.LpNorm;
    feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Targetting %f vertices in L%f norm", N, p);
    double n = 2.; // Space dimension
    double deg = (double)_options.polynomialDegree;

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

    metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, coord, exponentInIntegral, exponentForDeterminant);

    // Bound the eigenvalues according to hMin and hMax
    for(auto &pair : _metricTensorAtNodetags) {
      pair.second = pair.second.boundEigenvalues(_lambdaMin, _lambdaMax);
    }

    if(_options.debug)
      drawEllipsoids("metricsAfterScaling.pos", _metricTensorAtNodetags, nodeTags, coord, _options.plotSizeFactor, 20);
  }

  // Apply gradation
  // FIXME: gradation is computed for SMetric3 metric for now, so
  // need to transfer MetricTensors to SMetric3, then transfer back
  if(_options.enableGradation) {
    // applyGradation(nodeTags, coord);
    newGradation(nodeTags, coord, _metricTensorAtNodetags);
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