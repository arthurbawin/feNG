
#include "feMessage.h"
#include "feMetric.h"
#include "feMetricTools.h"
#include "feNumeric.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

extern int FE_VERBOSE;

#define MAX_VERTICES_BACKMESH 6

#if defined(HAVE_GMSH)
feStatus createMetricField(feMetric *&metricField,
                           feMetricOptions &options,
                           feMesh2DP1 *mesh,
                           std::vector<feNewRecovery*> recoveredFields)
{
  initializeGmsh();

  gmsh::open(options.backgroundMeshfile);
  gmsh::model::getCurrent(options.modelForMetric);

  if(metricField == nullptr) {
    metricField = new feMetric(options, mesh, recoveredFields);
  } else {
    feInfo("Pas nullptr");
    metricField->setOptions(options);
    metricField->setRecoveredFields(recoveredFields);
  }

  metricField->_backmeshOrder = mesh->getOrder();
  metricField->_nVerticesPerElmOnBackmesh = mesh->getNumVerticesPerElemForDimension(mesh->getDim());

  return FE_STATUS_OK;
}
#endif

feMetric::feMetric(const feMetricOptions &options,
                   feMesh2DP1 *mesh,
                   std::vector<feNewRecovery*> recoveredFields)
  : _options(options),
    _mesh(mesh),
    _recoveredFields(recoveredFields)
{
#if defined(HAVE_GMSH)
  // Min and max eigenvalues based on sizes
  _lambdaMax = 1. / (_options.hMin * _options.hMin);
  _lambdaMin = 1. / (_options.hMax * _options.hMax);

  // _nVerticesPerElm = mesh->getNumVerticesPerElemForDimension(mesh

  if(_recoveredFields.size() > 0)
    _options.polynomialDegree = _recoveredFields[0]->getDegreeSolution();

  // Create the Gmsh structures :
  // vertices tags and coordinates in the Gmsh model.
  // FIXME: Should eventually use internal structures only
  std::vector<double> parametricCoord;
  int dim = -1;
  int tag = -1;
  int includeBoundary = false;
  gmsh::model::mesh::getNodes(_nodeTags, _coord, parametricCoord, dim, tag, includeBoundary, false);

  // Create a nodeTags to feVertex->tag map (brute-force for now)
  const std::vector<Vertex> &meshVertices = mesh->getVertices();

  double tol = 1e-12;
  _v2n.clear();
  _n2v.clear();

  _sequentialTag2nodeTagVec.resize(meshVertices.size());

  for(size_t i = 0; i < _nodeTags.size(); i++) {
    const double x = _coord[3 * i + 0];
    const double y = _coord[3 * i + 1];
    for(const Vertex &v : meshVertices) {
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol) {
        _v2n[&v] = _nodeTags[i];
        _n2v[_nodeTags[i]] = &v;
        _nodeTag2sequentialTag[_nodeTags[i]] = mesh->_verticesMap[v.getTag()];
        _sequentialTag2nodeTag[mesh->_verticesMap[v.getTag()]] = _nodeTags[i];
        _sequentialTag2nodeTagVec[mesh->_verticesMap[v.getTag()]] = _nodeTags[i];
        break;
      }
    }
  }

  // Initialize the map/vector of metric tensors
  _metrics.resize(meshVertices.size());
  for(auto n : _nodeTags) {
    _metricTensorAtNodetags[n].setToIdentity();
  }
#endif
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

// Precompute all metric logarithms
void feMetric::computeLogMetrics()
{
  _logMetricTensorAtNodetags_eigen.resize(_metricTensorAtNodetags.size());
  for(auto &pair : _metricTensorAtNodetags) {
    _logMetricTensorAtNodetags[pair.first] = pair.second.log();
    _logMetricTensorAtNodetags_eigen[pair.first - 1] = convert2eigenMatrix(pair.second.log());
  }

  _logMetrics.resize(_metrics.size());
  for(size_t i = 0; i < _metrics.size(); ++i) {
    _logMetrics[i] = _metrics[i].log();
  }
}

// It's much faster to interpolate the determinant
// than interpolating the metric then taking the determinant.
double feMetric::computeIntegralOfDeterminant(const double exponentInIntegral)
{
  // tic();
  // double I = 0.;
  // #if defined(HAVE_OMP)
  // #pragma omp parallel
  // #endif
  // {
  //   double xsi[2];
  //   MetricTensor M_interpolated;
  //   int tags[MAX_VERTICES_BACKMESH];

  //   for(feMesh::iter_2D_Cnc it = _mesh->first2DCnc(); it != _mesh->last2DCnc(); ++it)
  //   {
  //     const feCncGeo *cnc = *(it);
  //     const int nElm = cnc->getNumElements();
  //     const int nQuad = cnc->getFeSpace()->getNumQuadPoints();
  //     const std::vector<double> &w = cnc->getFeSpace()->getQuadratureWeights();
  //     const std::vector<double> &xQ = cnc->getFeSpace()->getRQuadraturePoints();
  //     const std::vector<double> &yQ = cnc->getFeSpace()->getSQuadraturePoints();
  //     const std::vector<double> &J = cnc->getJacobians();

  //     #pragma omp for reduction(+ : I)
  //     for(int iElm = 0; iElm < nElm; ++iElm) {

  //       for(int i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
  //         tags[i] = cnc->getVertexConnectivity(iElm, i);
  //       }

  //       const MetricTensor &logM0 = _logMetrics[tags[0]];
  //       const MetricTensor &logM1 = _logMetrics[tags[1]];
  //       const MetricTensor &logM2 = _logMetrics[tags[2]];

  //       for(int k = 0; k < nQuad; ++k) {

  //         xsi[0] = xQ[k];
  //         xsi[1] = yQ[k];

  //         if(_nVerticesPerElmOnBackmesh == 3) {
  //           logEuclidianP1Interpolation(xsi, logM0, logM1, logM2, M_interpolated);
  //         }

  //         if(_nVerticesPerElmOnBackmesh == 6) {
  //           const MetricTensor &logM3 = _logMetrics[tags[3]];
  //           const MetricTensor &logM4 = _logMetrics[tags[4]];
  //           const MetricTensor &logM5 = _logMetrics[tags[5]];
  //           logEuclidianP2Interpolation(xsi, logM0, logM1, logM2, logM3, logM4, logM5,
  //                                       M_interpolated);
  //         }

  //         double interpolatedDet = M_interpolated.determinant();

  //         I += w[k] * J[iElm * nQuad + k] * pow(interpolatedDet, exponentInIntegral);
  //       }

  //     }

  //   }
  // }
  // feInfoCond(FE_VERBOSE > 0, "I1 = %+-1.16e in %f s", I, toc());

  // tic();
  double I = 0.;
  #if defined(HAVE_OMP)
  #pragma omp parallel
  #endif
  {
    double xsi[2];
    // MetricTensor M_interpolated;
    int tags[MAX_VERTICES_BACKMESH];

    for(feMesh::cncConstIter it = _mesh->first2DCnc(); it != _mesh->last2DCnc(); ++it)
    {
      const feCncGeo *cnc = *(it);
      const int nElm = cnc->getNumElements();
      const int nQuad = cnc->getFeSpace()->getNumQuadPoints();
      const std::vector<double> &w = cnc->getFeSpace()->getQuadratureWeights();
      const std::vector<double> &xQ = cnc->getFeSpace()->getRQuadraturePoints();
      const std::vector<double> &yQ = cnc->getFeSpace()->getSQuadraturePoints();
      const std::vector<double> &J = cnc->getJacobians();

      #if defined(HAVE_OMP)
      #pragma omp for reduction(+ : I)
      #endif
      for(int iElm = 0; iElm < nElm; ++iElm) {

        for(int i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
          tags[i] = cnc->getVertexConnectivity(iElm, i);
        }

        const double det0 = _metrics[tags[0]].determinant();
        const double det1 = _metrics[tags[1]].determinant();
        const double det2 = _metrics[tags[2]].determinant();

        for(int k = 0; k < nQuad; ++k) {
          xsi[0] = xQ[k];
          xsi[1] = yQ[k];
          double interpolatedDet = (1.-xsi[0]-xsi[1]) * det0 + xsi[0] * det1 + xsi[1] * det2;
          I += w[k] * J[iElm * nQuad + k] * pow(interpolatedDet, exponentInIntegral);
        }
      }
    }
  }
  // feInfoCond(FE_VERBOSE > 0, "I = %+-1.16e in %f s", I, toc());

  return I;
}

void feMetric::scaleMetricsByDeterminant(const double globalScalingFactor,
                                         const double exponentForDeterminant)
{
  for(auto &pair : _metricTensorAtNodetags) {
    const double det = pair.second.determinant();
    pair.second *= globalScalingFactor * pow(det, exponentForDeterminant);
  }
}

#if defined(HAVE_GMSH)
// Scale the metric field to fit N vertices in the final mesh. Interpolation is
// performed on the Gmsh substitute (gmshModel). Using gmsh api because it's easier.
// FIXME: parallelize
void feMetric::metricScalingFromGmshSubstitute(std::map<int, MetricTensor> &metrics,
                                               double exponentInIntegral,
                                               double exponentForDeterminant)
{

  // int nVerticesPerElement = _recoveredFields[0]->_cnc->getNumVerticesPerElem();

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

  // if(_nVerticesPerElmOnBackmesh == 3 && nVerticesPerElement == 6) {
  //   // Linear backmesh with a P2 feMesh : get P1 basis
  //   int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  //   gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents,
  //                                        basisFunctions, numOrientations);
  // }

  // Get the jacobians
  std::vector<double> jac, detJac, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, detJac, pts);

  // Compute integral of det^exponent
  double I = 0.0;
  size_t nQuad = weights.size();
  size_t nElm = elementTags[0].size();

  double N = (double)_options.nTargetVertices;
  double dim = 2.;

#if defined(HAVE_OMP)
#pragma omp parallel
#endif
  {
    MetricTensor M_interpolated;
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
  for(size_t i = 0; i < _nodeTags.size(); i++) {
    MetricTensor &M = metrics[_nodeTags[i]];
    double factor = pow(N / I, 2.0 / dim) * pow(M.determinant(), exponentForDeterminant);
    M *= factor;
  }
}
#endif

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

// Apply Alauzet's metric gradation
void feMetric::applyGradation()
{
  newGradation(_nodeTags, _coord, _options.gradation, _metricTensorAtNodetags);

  // FIXME: gradation is computed for SMetric3 metric for now, so
  // need to transfer MetricTensors to SMetric3, then transfer back
  // // Transfer from one map to the other to compute gradation
  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   _smetric3AtNodetags[nodeTags[i]] = convert2metric3(_metricTensorAtNodetags[nodeTags[i]]);
  // }

  // gradationMetriques(_options.gradation, 200, coord, _smetric3AtNodetags);

  // // Transfer back
  // for(size_t i = 0; i < nodeTags.size(); i++) {
  //   _metricTensorAtNodetags[nodeTags[i]](0, 0) = _smetric3AtNodetags[nodeTags[i]](0, 0);
  //   _metricTensorAtNodetags[nodeTags[i]](1, 0) = _smetric3AtNodetags[nodeTags[i]](1, 0);
  //   _metricTensorAtNodetags[nodeTags[i]](0, 1) = _smetric3AtNodetags[nodeTags[i]](0, 1);
  //   _metricTensorAtNodetags[nodeTags[i]](1, 1) = _smetric3AtNodetags[nodeTags[i]](1, 1);
  // }
}

#if defined(HAVE_GMSH)
void feMetric::writeMetricField()
{

  std::vector<std::size_t> P1ToP2;
  size_t numBGMVertices = _nodeTags.size();
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
      for(size_t j = 0; j < _nodeTags.size(); ++j) {
        double xp1 = coordP1[3 * i + 0];
        double yp1 = coordP1[3 * i + 1];
        double xp2 = _coord[3 * j + 0];
        double yp2 = _coord[3 * j + 1];
        if(sqrt((xp1 - xp2) * (xp1 - xp2) + (yp1 - yp2) * (yp1 - yp2)) < tol) {
          P1ToP2[i] = _nodeTags[j];
          break;
        }
      }
    }
  }

  std::string bgmToWrite;
  gmsh::model::getCurrent(bgmToWrite);

  // Create a view that contains the metric field
  _metricViewTag = gmsh::view::add(":metric");
  std::vector<std::vector<double> > metricData;

  std::vector<double> vMetric(9);

  for(size_t i = 0; i < numBGMVertices; i++) {

    // Get metric at vertex (if P1 background mesh) or at corresponding P2 vertex
    int tag = (_backmeshOrder > 1) ? P1ToP2[i] : _nodeTags[i];

    MetricTensor &M = _metricTensorAtNodetags[tag];

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.;

    metricData.push_back(vMetric);
  }

  if(_backmeshOrder > 1) {
    gmsh::view::addModelData(_metricViewTag, 0, bgmToWrite, "NodeData", nodeTagsP1, metricData);
  } else {
    gmsh::view::addModelData(_metricViewTag, 0, bgmToWrite, "NodeData", _nodeTags, metricData);
  }

  // MMG only takes .msh files of version 2.2
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  // Write mesh + metric
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);

  // Reset the initial model as active
  gmsh::model::setCurrent(initialMesh);
}
#endif

#if defined(HAVE_GMSH)
// Compute the metric tensor field from recovered derivatives,
// scale it to obtain roughly N vertices and write it as a Gmsh view
// to adapt with MMG
feStatus feMetric::computeMetrics()
{
  feInfoCond(FE_VERBOSE > 0, "");
  feInfoCond(FE_VERBOSE > 0, "RIEMANNIAN METRIC:");
  feInfoCond(FE_VERBOSE > 0, "\t\tReading mesh model with Gmsh");

  if(!_options.useAnalyticDerivatives && _recoveredFields.size() == 0)
  {
    return feErrorMsg(FE_STATUS_ERROR, "Provide either a solution recovery "
      "or analytical derivatives to compute a metric field.");
  }

  initializeGmsh();
  gmsh::open(_options.backgroundMeshfile);
  gmsh::model::getCurrent(_options.modelForMetric);

  feInfoCond(FE_VERBOSE >= 0, "\t\tComputing metrics on Gmsh model %s",
    _options.modelForMetric.data());

  double N = (double)_options.nTargetVertices;
  double p = _options.LpNorm;
  double n = 2.; // Space dimension
  double deg = (double)_options.polynomialDegree;

  // Compute the raw metric field (no scaling)
  switch(_options.method) {

    case adaptationMethod::ANALYTIC_METRIC:
      feCheckReturn(setAnalyticMetric(_nodeTags, _coord));
      _scaleMetricField = false;
      break;

    case adaptationMethod::ISO_P1: {
      bool isotropic = true;
      feCheckReturn(computeMetricsP1(_nodeTags, _coord, isotropic));
      break;
    }

    case adaptationMethod::ANISO_P1:
      feCheckReturn(computeMetricsP1(_nodeTags, _coord));
      break;

    case adaptationMethod::ISO_PN: {
      bool isotropic = true;
      feCheckReturn(computeMetricsPn(_nodeTags, _coord, isotropic));
      break;
    }

    case adaptationMethod::ANISO_P2:
      feCheckReturn(computeMetricsP2(_nodeTags, _coord));
      break;

    case adaptationMethod::ANISO_P1_REFERENCE_SPACE:
      {
        int maxIter = 10;
        for(int iter = 0; iter < maxIter; ++iter) {
          // Compute Q_triangle in the reference space, then Q = M^(1/2) * Q_triangle * M^(1/2)
          bool setToIdentity = (iter == 0);
          feCheckReturn(computeMetricsP1_referenceSpace(_nodeTags, _coord, setToIdentity));
          std::string name1 = "Q_iter" + std::to_string(iter) + ".pos";
          std::string name2 = "M_iter" + std::to_string(iter) + ".pos";
          drawEllipsoids(name1, _metricTensorAtNodetags, _nodeTags, _coord, _options.plotSizeFactor, 20);
          // Get M from Q
          // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
          // drawEllipsoids(name2, _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        }
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_final.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);

        // // Compute in the physical space to compare
        // s = computeMetricsP2(nodeTags, _coord);
        // drawEllipsoids("Q_physical.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_physical.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        exit(-1);
        break;
      }

    case adaptationMethod::ANISO_P2_REFERENCE_SPACE:
      {
        int maxIter = 10;
        for(int iter = 0; iter < maxIter; ++iter) {
          // Compute Q_triangle in the reference space, then Q = M^(1/2) * Q_triangle * M^(1/2)
          bool setToIdentity = (iter == 0);
          feCheckReturn(computeMetricsP2_referenceSpace(_nodeTags, _coord, setToIdentity));
          std::string name1 = "Q_iter" + std::to_string(iter) + ".pos";
          std::string name2 = "M_iter" + std::to_string(iter) + ".pos";
          // drawEllipsoids(name1, _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
          // Get M from Q
          // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
          // drawEllipsoids(name2, _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        }
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_final.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);

        // Compute in the physical space to compare
        feCheckReturn(computeMetricsP2(_nodeTags, _coord));
        // drawEllipsoids("Q_physical.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        // metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, _coord, exponentInIntegral, exponentForDeterminant);
        // drawEllipsoids("M_physical.pos", _metricTensorAtNodetags, nodeTags, _coord, _options.plotSizeFactor, 20);
        exit(-1);
        break;
      }

    case adaptationMethod::ANISO_PN:
        feCheckReturn(computeMetricsPn(_nodeTags, _coord));
      break;

    case adaptationMethod::GOAL_ORIENTED_ANISO_P1:
      feCheckReturn(computeMetricsGoalOrientedP1(_nodeTags, _coord));
      break;

    case adaptationMethod::CURVED_LS:
      feCheckReturn(computeMetricsCurvedLogSimplex(_nodeTags, _coord));
      break;

    case adaptationMethod::CURVED_LS_INDUCED_DIRECTIONS: {
      bool useInducedDirections = true;
      feCheckReturn(computeMetricsCurvedLogSimplex(_nodeTags, _coord, useInducedDirections));
      break;
    }

    case adaptationMethod::CURVED_ISOLINES:
      feCheckReturn(computeMetricsCurvedIsolines(_nodeTags, _coord));
      _scaleMetricField = _options.scaleIsolineMetric;
      break;

    case adaptationMethod::CURVED_REFERENCE_SPACE:
      feCheckReturn(computeMetricsCurvedReferenceSpace(_nodeTags, _coord));
      break;

    case adaptationMethod::CURVED_GRAPH_SURFACE:
      feCheckReturn(computeMetricsP2_forGraphSurface(_nodeTags, _coord));
      _scaleMetricField = false;
      break;
  }

  // Reconstruct the underlying manifold before scaling
  // reconstructManifoldFromMetric(_nodeTags, _coord, _metricTensorAtNodetags, _sequentialTag2nodeTagVec, _recoveredFields[0]);

  // Apply scaling and gradation only when not testing for metric interpolation convergence
  if(_options.debug)
    drawEllipsoids(_debugPrefix + "rawMetrics.pos", _metricTensorAtNodetags, _nodeTags, _coord, _options.plotSizeFactor, 20);

  if(_scaleMetricField)
  {
    double exponentInIntegral, exponentForDeterminant;

    if(_options.targetNorm == Norm::Lp) {
      exponentInIntegral = p * (deg + 1.0) / (2.* (p * (deg + 1.0) + n));
      exponentForDeterminant = -1. / (p * (deg + 1.) + n);
      feInfo("Exponents for scaling for L%2.1f norm = %+-1.4e - %+-1.4e", p, exponentForDeterminant, exponentInIntegral);
      feInfo("Target number of vertices = %+-1.4e", N);
    } else {
      exponentInIntegral = p * deg / (2.* (p * deg + n));
      exponentForDeterminant = -1. / (p * deg + n);
      feInfo("Exponents for scaling for H1 seminorm = %+-1.4e - %+-1.4e", exponentForDeterminant, exponentInIntegral);
    }

    // Scale the metric field
    feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Targetting %f vertices in L%f norm", N, p);
    metricScalingFromGmshSubstitute(_metricTensorAtNodetags, exponentInIntegral, exponentForDeterminant);

    // Bound the eigenvalues according to hMin and hMax
    for(auto &pair : _metricTensorAtNodetags) {
      pair.second = pair.second.boundEigenvalues(_lambdaMin, _lambdaMax);
    }

    if(_options.debug)
      drawEllipsoids(_debugPrefix + "metricsAfterScaling.pos", _metricTensorAtNodetags, _nodeTags, _coord, _options.plotSizeFactor, 20);
  }

  // Apply gradation
  if(_options.enableGradation) {
    newGradation(_nodeTags, _coord, _options.gradation, _metricTensorAtNodetags);
  }

  if(_options.debug) {
    drawEllipsoids(_debugPrefix + "metricsAfterGradation.pos",
      _metricTensorAtNodetags, _nodeTags, _coord, _options.plotSizeFactor, 20);
  }

  // Precompute all metric logarithms
  this->computeLogMetrics();

  // Write the metric field as a view in the Gmsh model
  writeMetricField();

  return FE_STATUS_OK;
}
#endif