set(feNG_SRC
  feBilinearForm.cpp
  feLinearSystem.cpp
  feMesh.cpp
  feNumber.cpp
  feSolution.cpp
  feSolutionContainer.cpp
  feSolver.cpp
  feSpace.cpp
  feSysElm.cpp
)

set(feNG_API
  feBilinearForm.h
  feCncGeo.h
  feFunction.h
  feLinearSystem.h
  feMesh.h
  feNG.h
  feNorm.h
  feNumber.h
  feQuadrature.h
  feSolution.h
  feSolutionContainer.h
  feSolver.h
  feSpace.h
  feSysElm.h
)

feNG_add_library(feNG "${feNG_SRC}" "${feNG_API}")
