set(feNG_SRC
  ellipseToolbox.cpp
  feBilinearForm.cpp
  feCncGeo.cpp
  feCompressedRowStorage.cpp
  # feCurvedAdapt.cpp
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
  feTriangle.cpp
)

if(ENABLE_EIGEN)
  # Add Gmsh files for STensor3
  list(APPEND feNG_SRC
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/fullMatrix.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/Numeric.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/STensor3.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/poly34.cpp)
endif()

set(feNG_API
  ellipseToolbox.h
  feBilinearForm.h
  feCncGeo.h
  feCompressedRowStorage.h
  # feCurvedAdapt.h
  # feEdge.h
  # feElement.h
  feExporter.h
  feFunction.h
  # feLine.h
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
  feTriangle.h
  feVertex.h
)

feNG_add_library(feNG "${feNG_SRC}" "${feNG_API}")