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
#include "feBilinearForm.h"
#include "feTimeIntegration.h"
#include "feNonLinearSolver.h"
#include "feLinearSystem.h"
#include "feEigenProblem.h"
#include "feExporter.h"
#include "feMetric.h"
#include "feRecovery.h"
#include "feNewRecovery.h"
#include "feNorm.h"
#include "feGraphics.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

#endif