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
#include "feNumeric.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

// Useful scalar and vector constants
namespace scalarConstant {
	feConstantFunction zero(0.);
  feConstantFunction one(1.);
  feConstantFunction minusOne(-1.);
}

namespace vectorConstant {
  feConstantVectorFunction zero({0., 0., 0.});
  feConstantVectorFunction one({1., 1., 1.});
  feConstantVectorFunction one_x({1., 0., 0.});
  feConstantVectorFunction one_y({0., 1., 0.});
  feConstantVectorFunction one_z({0., 0., 1.});
  feConstantVectorFunction minusOne_x({-1., 0., 0.});
  feConstantVectorFunction minusOne_y({0., -1., 0.});
  feConstantVectorFunction minusOne_z({0., 0., -1.});
}

#endif