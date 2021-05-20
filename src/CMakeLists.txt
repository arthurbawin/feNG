set(feNG_SRC
  feBilinearForm.cpp
  feCompressedRowStorage.cpp
  feExporter.cpp
  feLinearSystem.cpp
  feLinearSystemMklPardiso.cpp
  feLinearSystemPETSc.cpp
  feMesh.cpp
  feMeshRead.cpp
  feNorm.cpp
  feNumber.cpp
  feSolution.cpp
  feQuadrature2.cpp
  feSolutionContainer.cpp
  feSolver.cpp
  feSpace.cpp
  feSpaceTriangle.cpp
  feSysElm.cpp
)

set(feNG_API
  feBilinearForm.h
  feCncGeo.h
  feCompressedRowStorage.h
  feEdge.h
  feElement.h
  feExporter.h
  feFunction.h
  feLine.h
  feLinearSystem.h
  feLinearSystemMklPardiso.h
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
