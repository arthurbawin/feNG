set(feNG_SRC
  feBilinearForm.cpp
  feMesh.cpp
  feNumber.cpp
  feSolution.cpp
  feSpace.cpp
  feSysElm.cpp
  feTimeIntegrator.cpp
)

set(feNG_API
  feBilinearForm.h
  feCncGeo.h
  feFunction.h
  feLinearSystem.h
  feMesh.h
  feNG.h
  feNumber.h
  feQuadrature.h
  feSolution.h
  feSpace.h
  feSysElm.h
  feTimeIntegrator.h
)

feNG_add_library(feNG "${feNG_SRC}" "${feNG_API}")
