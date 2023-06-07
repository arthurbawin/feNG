
#include "ellipseToolbox.h"
#include "fullMatrix.h"
#include "feMessage.h"
#include "feMetric.h"
#include "feMetricTools.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif
#if defined(HAVE_SOPLEX)
#include "soplex.h"
using namespace soplex;
#endif

extern int FE_VERBOSE;

// #define ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE

feMetric::feMetric(feRecovery *recovery, feMetricOptions metricOptions)
  : _recovery(recovery), _options(metricOptions)
{
  _options.polynomialDegree = _recovery->getDegreeSolution();
}

feMetric::feMetric(feNewRecovery *recovery, feMetricOptions metricOptions)
  : _newRecovery(recovery), _options(metricOptions)
{
  _options.polynomialDegree = _newRecovery->getDegreeSolution();
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
  if(_newRecovery != nullptr)
    mesh = _newRecovery->_mesh;
  else if(_recovery != nullptr)
    mesh = _recovery->_mesh;
  else
    return feErrorMsg(FE_STATUS_ERROR, "No valid recovery");

  std::vector<Vertex> &meshVertices = mesh->getVertices();

  double tol = 1e-5;
  _v2n.clear();
  _n2v.clear();
  for(size_t i = 0; i < nodeTags.size(); i++) {
    const double x = coord[3 * i + 0];
    const double y = coord[3 * i + 1];
    for(auto &v : meshVertices) {
      if(fabs(x - v.x()) < tol && fabs(y - v.y()) < tol) {
        _v2n[&v] = nodeTags[i];
        _n2v[nodeTags[i]] = &v;
        _nodeTag2sequentialTag[nodeTags[i]] = mesh->_verticesMap[v.getTag()];
        _sequentialTag2nodeTag[mesh->_verticesMap[v.getTag()]] = nodeTags[i];
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
                                               double exponentInIntegral,
                                               double exponentForDeterminant)
{
#if defined(HAVE_GMSH)

	int nVerticesPerElement = _newRecovery->_cnc->getNumVerticesPerElem();

  // Get the mesh elements
  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t> > elementTags;
  std::vector<std::vector<std::size_t> > elemNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);

  if(elementTypes.size() > 1){
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
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents, basisFunctions,
                                       numOrientations);

  if(_nVerticesPerElmOnBackmesh == 3 && nVerticesPerElement == 6) {
  	// Linear backmesh with a P2 feMesh : get P1 basis
  	int triP1 = gmsh::model::mesh::getElementType("Triangle", 1);
  	gmsh::model::mesh::getBasisFunctions(triP1, localCoord, "Lagrange", numComponents, basisFunctions,
                                       numOrientations);
  }

  // Get the jacobians
  std::vector<double> jac, det, pts;
  gmsh::model::mesh::getJacobians(elementType, localCoord, jac, det, pts);

  // Compute integral of det^exponent
  double I = 0.0;
  int nQuad = weights.size();
  size_t nElm = elementTags[0].size();

  double N = (double) _options.nTargetVertices;
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
    #pragma omp for reduction(+:I)
    #endif
    for(size_t iElm = 0; iElm < nElm; iElm++) {

    	for(size_t i = 0; i < _nVerticesPerElmOnBackmesh; ++i) {
    		tags[i] = elemNodeTags[0][_nVerticesPerElmOnBackmesh * iElm + i];
    	}

      for(size_t i = 0; i < nQuad; i++) {

      	xsi[0] = localCoord[3 * i + 0];
      	xsi[1] = localCoord[3 * i + 1];

      	if(_nVerticesPerElmOnBackmesh == 3){
      		MetricTensor &M0 = metrics.at(tags[0]);
      		MetricTensor &M1 = metrics.at(tags[1]);
      		MetricTensor &M2 = metrics.at(tags[2]);
      		logEuclidianP1Interpolation(xsi, M0, M1, M2, M_interpolated);
      	}

      	if(_nVerticesPerElmOnBackmesh == 6){
      		MetricTensor &M0 = metrics.at(tags[0]);
      		MetricTensor &M1 = metrics.at(tags[1]);
      		MetricTensor &M2 = metrics.at(tags[2]);
      		MetricTensor &M3 = metrics.at(tags[3]);
      		MetricTensor &M4 = metrics.at(tags[4]);
      		MetricTensor &M5 = metrics.at(tags[5]);
      		logEuclidianP2Interpolation(xsi, M0, M1, M2, M3, M4, M5, M_interpolated);
      	}

      	double interpolatedDet = M_interpolated.determinant();

      	if(interpolatedDet <= 0. ) {
          feInfo("Negative metric determinant:");
      		M_interpolated.print();
      		exit(-1);
      	}

        I += weights[i] * det[iElm * nQuad + i] * pow( interpolatedDet, exponentInIntegral);
      }
    }
  }

  feInfo("Scaling : Integral of determinant = %+-1.16e", I);

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

inline double matfun(double x, double y){
	return 1. + x*x*x + y*y;
}

MetricTensor analyticMetric(double x, double y)
{
	MetricTensor res(1.0);
	SVector3 EX(1., 0., 0.);
	SVector3 EY(0., 1., 0.);

	SVector3 e1 = EX *   cos(M_PI*x)  + EY * sin(M_PI*x);
	SVector3 e2 = EX * (-sin(M_PI*x)) + EY * cos(M_PI*x);
	SVector3 e3(0., 0., 1.);
	double h1 = 1 + x*y;
	double h2 = 5 + x*y*y;
	SMetric3 tmp(1./(h1*h1), 1./(h2*h2), 1., e1, e2, e3);
	res(0,0) = tmp(0,0);
	res(0,1) = tmp(0,1);
	res(1,0) = tmp(1,0);
	res(1,1) = tmp(1,1);
	return res;
}

void feMetric::interpolationTest(const std::vector<size_t> &nodeTags, std::vector<double> &coord)
{
#if defined(HAVE_GMSH)

	int nVerticesPerElement = _newRecovery->_cnc->getNumVerticesPerElem();

  // Get quadrature rule and interpolation functions on the Gmsh substitute
  int elementOrder = getGeometricInterpolantDegree(_newRecovery->_cnc->getInterpolant());
  feInfo("Mesh elements order: %d", elementOrder);
  int elementType = gmsh::model::mesh::getElementType("Triangle", elementOrder);

  int numComponents, numOrientations;
  std::vector<double> localCoord;
  std::vector<double> weights;
  std::vector<double> basisFunctions;
  gmsh::model::mesh::getIntegrationPoints(elementType, "Gauss12", localCoord, weights);
  gmsh::model::mesh::getBasisFunctions(elementType, localCoord, "Lagrange", numComponents, basisFunctions,
                                       numOrientations);

  // Get the mesh elements
  // std::vector<int> elementTypes;
  // std::vector<std::vector<std::size_t> > elementTags;
  // std::vector<std::vector<std::size_t> > elemNodeTags;
  // gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, 2);
  std::vector<std::size_t> elementTags;
  std::vector<std::size_t> elemNodeTags;
  gmsh::model::mesh::getElementsByType(elementType, elementTags, elemNodeTags);

  if(elementTags.empty()) {
    feErrorMsg(FE_STATUS_ERROR, "There are no P2 elements in the Gmsh model. It's probably because the Gmsh model is still P1 for MMG.");
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

  for(size_t iElm = 0; iElm < elementTags.size(); iElm++)
  {
  	int v0 = elemNodeTags[nVerticesPerElement * iElm + 0] - 1;
  	int v1 = elemNodeTags[nVerticesPerElement * iElm + 1] - 1;
  	int v2 = elemNodeTags[nVerticesPerElement * iElm + 2] - 1;
  	int v3 = elemNodeTags[nVerticesPerElement * iElm + 3] - 1;
  	int v4 = elemNodeTags[nVerticesPerElement * iElm + 4] - 1;
  	int v5 = elemNodeTags[nVerticesPerElement * iElm + 5] - 1;

  	double x0 = coord[3 * v0];
  	double y0 = coord[3 * v0+1];
  	double x1 = coord[3 * v1];
  	double y1 = coord[3 * v1+1];
  	double x2 = coord[3 * v2];
  	double y2 = coord[3 * v2+1];
  	double x3 = coord[3 * v3];
  	double y3 = coord[3 * v3+1];
  	double x4 = coord[3 * v4];
  	double y4 = coord[3 * v4+1];
  	double x5 = coord[3 * v5];
  	double y5 = coord[3 * v5+1];

  	if(iElm == 0){
  		fprintf(myfile, "SP(%f,%f,0.){0.};\n", x0, y0);
  		fprintf(myfile, "SP(%f,%f,0.){1.};\n", x1, y1);
  		fprintf(myfile, "SP(%f,%f,0.){2.};\n", x2, y2);
  		fprintf(myfile, "SP(%f,%f,0.){3.};\n", x3, y3);
  		fprintf(myfile, "SP(%f,%f,0.){4.};\n", x4, y4);
  		fprintf(myfile, "SP(%f,%f,0.){5.};\n", x5, y5);
  	}

  	double X[6] = {x0, x1, x2, x3, x4, x5};
  	double Y[6] = {y0, y1, y2, y3, y4, y5};
  	double U[6] = {matfun(x0,y0), matfun(x1,y1), matfun(x2,y2), matfun(x3,y3), matfun(x4,y4), matfun(x5,y5)};

  	// MetricTensor M0(1.0); M0(0,0) = matfun(x0,y0); M0(1,1) = matfun(x0,y0);
  	// MetricTensor M1(1.0); M1(0,0) = matfun(x1,y1); M1(1,1) = matfun(x1,y1);
  	// MetricTensor M2(1.0); M2(0,0) = matfun(x2,y2); M2(1,1) = matfun(x2,y2);
  	// MetricTensor M3(1.0); M3(0,0) = matfun(x3,y3); M3(1,1) = matfun(x3,y3);
  	// MetricTensor M4(1.0); M4(0,0) = matfun(x4,y4); M4(1,1) = matfun(x4,y4);
  	// MetricTensor M5(1.0); M5(0,0) = matfun(x5,y5); M5(1,1) = matfun(x5,y5);

  	MetricTensor M0 = analyticMetric(x0,y0);
  	MetricTensor M1 = analyticMetric(x1,y1);
  	MetricTensor M2 = analyticMetric(x2,y2);
  	MetricTensor M3 = analyticMetric(x3,y3);
  	MetricTensor M4 = analyticMetric(x4,y4);
  	MetricTensor M5 = analyticMetric(x5,y5);

  	if(iElm == 0){
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
    	double sum  = 0.;
    	for(int ii = 0; ii < 6; ++ii){
    		sum   += basisFunctions[nVerticesPerElement * i + ii];
    		xphys += basisFunctions[nVerticesPerElement * i + ii] * X[ii];
    		yphys += basisFunctions[nVerticesPerElement * i + ii] * Y[ii];
    		uh    += basisFunctions[nVerticesPerElement * i + ii] * U[ii];
    	}

    	// MetricTensor Mkexact(1.0);
    	// Mkexact(0,0) = matfun(xphys,yphys);
    	// Mkexact(1,1) = matfun(xphys,yphys);
    	MetricTensor Mkexact = analyticMetric(xphys,yphys);

    	if(iElm == 0){
    		double POS_QUAD[2] = {xphys, yphys};
    		drawSingleEllipse(myfile2, POS_QUAD, Mkexact, scale, 100);
    		drawSingleEllipse(myfile3, POS_QUAD, Mk, scale, 100);
    	}

    	double uexact = matfun(xphys,yphys);

    	double frobError = matNorm2(Mk, Mkexact);

    	// Erreur sur l'interpolation des metriques
    	error += frobError*frobError * weights[i] * det[iElm * nQuad + i];

    	// Erreur (test de verification) sur l'interpolation d'un champ
    	// error += (uh-uexact)*(uh-uexact) * weights[i] * det[iElm * nQuad + i];

      area += weights[i] * det[iElm * nQuad + i];
    }
  }

  fprintf(myfile, "};\n"); fclose(myfile);
  fprintf(myfile2, "};\n"); fclose(myfile2);
  fprintf(myfile3, "};\n"); fclose(myfile3);
  fprintf(myfile4, "};\n"); fclose(myfile4);

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

    MetricTensor M = analyticMetric(x[0],x[1]);

    drawSingleEllipse(file, x, M, 500., 30);
  }

  fprintf(file, "};\n"); fclose(file);
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

// template <class MetricType> Eigen::Matrix2d convert2eigenMatrix(const MetricType &other)
// {
//   Eigen::Matrix2d M = Eigen::Matrix2d::Identity();
//   M(0, 0) = other(0, 0);
//   M(1, 0) = other(1, 0);
//   M(0, 1) = other(0, 1);
//   M(1, 1) = other(1, 1);
//   return M;
// }

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
  // Create a view that contains the metric field
  _metricViewTag = gmsh::view::add(":metric");
  int recoveryViewTag = gmsh::view::add(":recovery");
  std::vector<std::vector<double> > metricData;
  std::vector<std::vector<double> > recoveryData;

  double x[2];
  std::vector<double> vMetric(9);
  std::vector<double> vRecovery(1);
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    MetricTensor &M = _metricTensorAtNodetags[nodeTags[i]];

    // For interpolation in curved mesh (still requires Eigen matrices)
    // REMOVED TO USE LOGEUCLIDIAN INTERPOLATION?
    // _metricsOnGmshModel_eigen[nodeTags[i]] = convert2eigenMatrix(M);

    vMetric[0] = M(0, 0);
    vMetric[1] = M(0, 1);
    vMetric[2] = 0;

    vMetric[3] = M(0, 1);
    vMetric[4] = M(1, 1);
    vMetric[5] = 0;

    vMetric[6] = 0;
    vMetric[7] = 0;
    vMetric[8] = 1.0;
    vRecovery[0] = _newRecovery->evaluateRecovery(PPR::RECOVERY, 0, x);

    metricData.push_back(vMetric);
    recoveryData.push_back(vRecovery);
  }

  gmsh::view::addModelData(_metricViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           metricData);
  gmsh::view::addModelData(recoveryViewTag, 0, _options.modelForMetric, "NodeData", nodeTags,
                           recoveryData);

  // MMG only take .msh files of version 2.2
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  // Write mesh + metric and mesh + recovery
  gmsh::view::write(_metricViewTag, _options.mmgInputMeshfile);
  gmsh::view::write(recoveryViewTag, _options.recoveryName);
#endif
}

// Compute scaled optimal metric field for P1 elements
// according to Alauzet & Loseille
feStatus feMetric::computeMetricsP1(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
  // Min and max eigenvalues based on sizes
  double lambdaMax = 1. / (_options.hMin * _options.hMin);
  double lambdaMin = 1. / (_options.hMax * _options.hMax);

  double x[2], fxx, fxy, fyx, fyy;
  MetricTensor H;

  // Compute bounded absolute value of hessian at vertices (at nodetags)
  for(size_t i = 0; i < nodeTags.size(); i++) {
    x[0] = coord[3 * i + 0];
    x[1] = coord[3 * i + 1];

    // Evaluate 2nd order derivatives
    fxx = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 2, x);
    fxy = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 3, x);
    fyx = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 4, x);
    fyy = _newRecovery->evaluateRecovery(PPR::DERIVATIVE, 5, x);

    H(0, 0) = fxx;
    H(0, 1) = (fxy + fyx) / 2.;
    H(1, 0) = (fxy + fyx) / 2.;
    H(1, 1) = fyy;
    _metricTensorAtNodetags[nodeTags[i]] = H.boundEigenvaluesOfAbs(lambdaMin, lambdaMax);
  }

  return FE_STATUS_OK;
}

#if defined(HAVE_SOPLEX)
void setUpLinearProblem(linearProblem &myLP, feMetricOptions &options, int nTheta, bool reset = false)
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
  myLP.lprow = LPRow(3);
  myLP.lprow.setRhs(infty);
  myLP.lprowset = LPRowSet(numConstraints, 3);

  for(int i = 0; i < numConstraints; ++i) {
    DSVector row(3);
    row.add(0, 1.0);
    row.add(1, 1.0);
    row.add(2, 1.0);
    myLP.problem.addRowReal(LPRow(1.0, row, infty));
  }
}
#endif

// Compute scaled optimal metric field for Pn elements
// according to Coulaud & Loseille using the log-simplex method
feStatus feMetric::computeMetricsPn(std::vector<std::size_t> &nodeTags, std::vector<double> &coord)
{
#if defined(ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE)
  // Just assign analytic metric and exit
  // To test convergence of metric interpolation
  for(size_t i = 0; i < numVertices; i++) {
    double X = coord[3 * i + 0];
    double Y = coord[3 * i + 1];
    _metricTensorAtNodetags[nodeTags[i]] = analyticMetric(X, Y);
  }
  return FE_STATUS_OK;
#endif

#if defined(HAVE_SOPLEX)

  std::map<int, std::vector<double> > &errorCoeffAtVertices = _newRecovery->getErrorCoefficients();

  bool OK = true;

  size_t numVertices = nodeTags.size();
  int cnter = 0, maxThreads = omp_get_max_threads();

  #if defined(HAVE_OMP)
  #pragma omp parallel shared(cnter, OK)
  #endif
  {
      // Min and max eigenvalues based on sizes
    double lambdaMax = 1. / (_options.hMin * _options.hMin);
    double lambdaMin = 1. / (_options.hMax * _options.hMax);
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
    #pragma omp for schedule(dynamic,1)
    #endif
	  for(size_t i = 0; i < numVertices; i++) {

      if(!OK) continue;

	    x[0] = coord[3 * i + 0];
	    x[1] = coord[3 * i + 1];

	    // Get the coefficients of the homogeneous error polynomial at vertex
	    std::vector<double> &errorCoeff = errorCoeffAtVertices[_nodeTag2sequentialTag[nodeTags[i]]];

	    // Mixed derivatives coefficient
	    if(_options.polynomialDegree == 1){
	    	errorCoeff[1] /= 2.;
	    } else if(_options.polynomialDegree == 2){
	    	errorCoeff[1] /= 3.;
	    	errorCoeff[2] /= 3.;
	    } else {
	    	// return feErrorMsg(FE_STATUS_ERROR, "Modify mixed derivative error coefficients");
	    }

		  bool res = computeMetricLogSimplexStraight(x, errorCoeff, _options.polynomialDegree, nTheta,
		                                               maxIter, tol, Q, numIter, myLP);

      // int retry = 0, maxTry = 1;
      // while(!res && retry < maxTry) {
      //   nTheta *= 2;
      //   feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Could not compute metric. Trying again with %d constraints.", nTheta);
      //   setUpLinearProblem(myLP, _options, nTheta);
      //   res = computeMetricLogSimplexStraight(x, errorCoeff, _options.polynomialDegree, nTheta,
      //                                              maxIter, tol, Q, numIter, myLP);
      // }

      // Restore initial number of constraints
      // nTheta = _options.logSimplexOptions.nThetaPerQuadrant;

	    if(res) {
	      feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, 
	      			 "Computed metric in %2d iterations  - vertex %6d/%6d",
	             numIter, ++cnter, numVertices);

	      // Bound the eigenvalues of Q
	      #if defined(HAVE_OMP)
			  #pragma omp critical
			  #endif
	      _metricTensorAtNodetags[nodeTags[i]] = Q.boundEigenvaluesOfAbs(lambdaMin, lambdaMax);

	    } else {
	    	OK = false;
	      // return feErrorMsg(FE_STATUS_ERROR, "Could not compute a metric at (%+-1.5e - %+-1.5e) (vertex %d/%d)",
	      //           x[0], x[1], i, nodeTags.size());
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

// Compute the metric tensor field from recovered derivatives,
// scale it to obtain roughly N vertices and write it as a Gmsh view
// to adapt with MMG
feStatus feMetric::computeMetrics()
{
#if defined(HAVE_GMSH)

  bool debug = false;

  if(!_options.isGmshModelReady) {
    return feErrorMsg(FE_STATUS_ERROR, "No Gmsh model available to compute a metric field.\n"
                                       "Create a Gmsh mesh first.\n");
  }

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
    case adaptationMethod::ANISO_P1:
      s = computeMetricsP1(nodeTags, coord);
      break;

    case adaptationMethod::ANISO_PN:
      s = computeMetricsPn(nodeTags, coord);
      break;

      // case adaptationMethod::CURVED_LS:
      //   return computeMetricsLogSimplex();

      // case adaptationMethod::CURVED_EXTREME_SIZES:
      //   return computeMetricsExtremeSizesOnly();
  }

  if(s != FE_STATUS_OK) return s;

  #if !defined(ONLY_TEST_METRIC_INTERPOLATION_CONVERGENCE)
  // Apply scaling and gradation only when not testing for metric interpolation convergence
  if(debug) drawEllipsoids("rawMetrics.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);

  // Scale the metric field
  double N = (double) _options.nTargetVertices;
  double p = _options.LpNorm;
  feInfoCond(FE_VERBOSE >= VERBOSE_MODERATE, "Targetting %f vertices in L%f norm", N, p);
  double n = 2.; // Space dimension
  double deg = (double) _options.polynomialDegree;

  double exponentInIntegral, exponentForDeterminant;
  switch(_options.method) {
    case adaptationMethod::ANISO_P1:
      exponentInIntegral = p / (2. * p + n);
      exponentForDeterminant = -1. / (2. * p + n);
      break;
    case adaptationMethod::ANISO_PN:
      // exponentInIntegral = p * (deg + 1.0) / (p * (deg + 1.0) + n);
      exponentInIntegral = p * ((deg + 1.0)/2.) / (p * (deg + 1.0) + n);
      exponentForDeterminant = -1. / (p * (deg + 1.) + n);
      break;
     default:
     	return feErrorMsg(FE_STATUS_ERROR, "No exponents provided to scale the metric field");
  }

  metricScalingFromGmshSubstitute(_metricTensorAtNodetags, nodeTags, exponentInIntegral,
                                  exponentForDeterminant);

  if(debug)
    drawEllipsoids("metricsAfterScaling.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);

  // Apply gradation
  // FIXME: gradation is computed for SMetric3 metric for now, so
  // need to transfer MetricTensors to SMetric3, then transfer back
  if(_options.enableGradation) {
    applyGradation(nodeTags, coord);
  }

  if(debug)
    drawEllipsoids("metricsAfterGradation.pos", _metricTensorAtNodetags, nodeTags, coord, 100, 30);
  #endif

  // Write the metric field as a view in the Gmsh model
  writeMetricField(nodeTags, coord);

#endif
  return FE_STATUS_OK;
}