#ifndef _FEAPI_
#define _FEAPI_

#include "feNG.h"
#include "feMessage.h"
#include "feFunction.h"
#include "feOptionsParser.h"

#include "feCncGeo.h"
#include "feMesh.h"
#include "feQuadrature.h"
#include "feSpace.h"
#include "feNumber.h"
#include "feSolution.h"
#include "feSysElm.h"
#include "feSysElm_Ext.h"
#include "feBilinearForm.h"
#include "feSolver.h"
#include "feLinearSystem.h"
#include "feEigenProblem.h"
#include "feExporter.h"
#include "feRecovery.h"
#include "feNewRecovery.h"
#include "feMetric.h"
#include "feCurvedAdapt.h"
#include "feAdaptMesh.h"
#include "feNorm.h"
#include "feGraphics.h"

#if defined(HAVE_OMP)
#include "omp.h"
#endif

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

#endif