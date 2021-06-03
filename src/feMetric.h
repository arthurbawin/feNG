#ifndef _FEMETRIC_
#define _FEMETRIC_

#include "feRecovery.h"

#include "STensor3.h"

typedef struct feMetricOptions{
  int nTargetVertices;
  double hMin;
  double hMax;
  int computationMethod;
  double LpNorm;
} feMetricOptions;

class feMetric{
protected:
	feRecovery *_recovery;
  feMetricOptions _options;

  std::map<int, SMetric3> _metrics;
public:
  feMetric(feRecovery *recovery, feMetricOptions metricOptions);
  ~feMetric(){
  }

  void metricScaling();
  void metriqueSimplexe2D(int nPhi, std::vector<double> phi, std::vector<double> erreur, double &A, double &B, double &C, int max_iter);

  // Intersecter avec un autre feMetric
  void writeSizeFieldSol(std::string solFileName);
  void writeSizeFieldGmsh(std::string meshName, std::string metricMeshName);

};

#endif