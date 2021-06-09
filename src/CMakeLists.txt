set(feNG_SRC
  feBilinearForm.cpp
  feCncGeo.cpp
  feCompressedRowStorage.cpp
  feExporter.cpp
  feLinearSystem.cpp
  feLinearSystemMklPardiso.cpp
  feLinearSystemPETSc.cpp
  feMesh.cpp
  feMeshRead.cpp
  feMetric.cpp
  feNorm.cpp
  feNumber.cpp
  feRecovery.cpp
  feQuadrature.cpp
  feSimplex.cpp
  feSolution.cpp
  feSolutionContainer.cpp
  feSolver.cpp
  feSpace.cpp
  feSpaceTriangle.cpp
  feSysElm.cpp
)

if(ENABLE_EIGEN)
  # Add Gmsh files for STensor3
  list(APPEND feNG_SRC
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/fullMatrix.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/Numeric.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/STensor3.cpp)
endif()

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
  feMetric.h
  feNG.h
  feNorm.h
  feNumber.h
  feQuadrature.h
  feRecovery.h
  feSimplex.h
  feSolution.h
  feSolutionContainer.h
  feSolver.h
  feSpace.h
  feSysElm.h
  feVertex.h
)

feNG_add_library(feNG "${feNG_SRC}" "${feNG_API}")

# if(ENABLE_PETSC)
#   list(APPEND feNG_SRC feLinearSystemPETSc.cpp)
#   list(APPEND feNG_API feLinearSystemPETSc.h)
# endif()
# if(ENABLE_MKL)
#   list(APPEND feNG_SRC feLinearSystemPardiso.cpp)
#   list(APPEND feNG_API feLinearSystemPardiso.h)
# endif()


