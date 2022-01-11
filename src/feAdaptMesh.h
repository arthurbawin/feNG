#ifndef _FEADAPTEMESH_
#define _FEADAPTEMESH_

#include "feMetric.h"

/* The meteaNumber, solution, interpolation space and recovery are needed to
define the error callback, the objective function to minimize. */
void createCurvedMesh(feFunction *solExact, feMetaNumber *metaNumber, feSolution *sol, feSpace *intSpace, feRecovery *recovery, feMetric *metric, feMetricOptions metricOptions);

#endif