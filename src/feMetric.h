#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"

#include "STensor3.h"

typedef struct feMetricOptions {
  bool debug = true;

  int computationMethod;
  int polynomialDegree;
  int nTargetVertices;
  double eTargetError;
  double hMin;
  double hMax;
  double LpNorm;
  int nPhi;

  // The name of the Gmsh model on which the metric field will be computed
  std::string gmshModel;
  bool isGmshModelReady = false;
  // Geometric data to create the model with Gmsh api
  double modelSizeX;
  double modelSizeY;
  int directionFieldFromDerivativesOfOrder;

  bool (*inside)(double *);
} feMetricOptions;

class feMetric {
protected:
  feRecovery *_recovery;
  feMetricOptions _options;

  std::map<int, SMetric3> _metrics;

  // The tag of the gmsh view in which the metric field is stored (if using Gmsh)
  int _metricViewTag = -1;

public:
  feMetric(feRecovery *recovery, feMetricOptions metricOptions);
  ~feMetric() {}

  void computeMetrics();
  void computeMetricsHechtKuate();
  void computeMetricsWithDirectionField();

  int getMetricViewTag(){ return _metricViewTag; }

  void metricScaling();
  // void metriqueSimplexe2D(int nPhi, std::vector<double> phi, std::vector<double> erreur, double &A,
  //                         double &B, double &C, int max_iter);

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol2D(std::string solFileName);
  void writeSizeFieldSol3D(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);
  void drawEllipsoids(std::string posFile);
};

#endif