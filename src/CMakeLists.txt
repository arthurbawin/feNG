set(feNG_SRC
  feBilinearForm.cpp
  feCompressedRowStorage.cpp
  feLinearSystem.cpp
  feLinearSystemPardiso.cpp
  feLinearSystemPETSc.cpp
  feMesh.cpp
  feNorm.cpp
  feNumber.cpp
  feSolution.cpp
  feQuadrature2.cpp
  feSolutionContainer.cpp
  feSolver.cpp
  feSpace.cpp
  feSysElm.cpp
)

set(feNG_API
  feBilinearForm.h
  feCncGeo.h
  feCompressedRowStorage.h
  feFunction.h
  feLinearSystem.h
  feLinearSystemPardiso.h
  feLinearSystemPETSc.h
  feMesh.h
  feNG.h
  feNorm.h
  feNumber.h
  feQuadrature.h
  feQuadrature2.h
  feSolution.h
  feSolutionContainer.h
  feSolver.h
  feSpace.h
  feSysElm.h
  
)

# if(ENABLE_PETSC)
#   list(APPEND feNG_SRC feLinearSystemPETSc.cpp)
#   list(APPEND feNG_API feLinearSystemPETSc.h)
# endif()
# if(ENABLE_MKL)
#   list(APPEND feNG_SRC feLinearSystemPardiso.cpp)
#   list(APPEND feNG_API feLinearSystemPardiso.h)
# endif()

feNG_add_library(feNG "${feNG_SRC}" "${feNG_API}")
