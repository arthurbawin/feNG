#ifndef _FEAPI_
#define _FEAPI_

#include "feBilinearForm.h"
#include "feCncGeo.h"
#include "feConstants.h"
#include "feEigenProblem.h"
#include "feExporter.h"
#include "feFunction.h"
#include "feGraphics.h"
#include "feLinearSystem.h"
#include "feMesh.h"
#include "feMessage.h"
#include "feMetric.h"
#include "feNG.h"
#include "feNewRecovery.h"
#include "feNonLinearSolver.h"
#include "feNorm.h"
#include "feNumber.h"
#include "feNumeric.h"
#include "feOptionsParser.h"
#include "feQuadrature.h"
#include "feRecovery.h"
#include "feSolution.h"
#include "feSpace.h"
#include "feSysElm.h"
#include "feTimeIntegration.h"

#if defined(HAVE_GMSH)
#  include "gmsh.h"
#endif

#endif