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
#include "feSolver.h"
#include "feLinearSystem.h"
#include "feLinearSystemPETSc.h"
#include "feExporter.h"
#include "feRecovery.h"
#include "feMetric.h"
#include "feCurvedAdapt.h"
#include "feAdaptMesh.h"
#include "fePostProc.h"
#include "feNorm.h"

#if defined(HAVE_GMSH)
#include "gmsh.h"
#endif

#endif