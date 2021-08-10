#ifndef _FECURVEDADAPT_
#define _FECURVEDADAPT_

#include "feRecovery.h"
#include "feMetric.h"
// #include "feGridFunction.h"

#include "rtree.h"

class feCurvedAdapt {
protected:
  std::vector<feRecovery *> &_rec;
  // std::vector<feGridFunction*> _derivatives;

  RTree<int, double, 3> _rtree;

private:
  // void computeMetricP2(int node, double lMin, double lMax, double eps, double &g00, double &g01,
  // double &g11, FILE *F, double &C, double &S);
  void computeMetricP2(double x, double y, double lMin, double lMax, double eps, double &g00,
                       double &g01, double &g11, FILE *F, double &C, double &S);

public:
  feCurvedAdapt(feMesh *mesh, std::vector<feRecovery *> &recovery, feMetricOptions metricOptions, std::string metricMeshName = "", std::string nextMeshName = "", int useAnalytical = 0);
  ~feCurvedAdapt();
};

#endif