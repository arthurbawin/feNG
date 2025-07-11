set(feNG_SRC
    CHNS_Solver.cpp
    DiffusionSolver.cpp
    ellipseToolbox.cpp
    feAdaptMesh.cpp
    feBilinearForm.cpp
    feBoundaryConditions.cpp
    feCncGeo.cpp
    feColoring.cpp
    feCompressedRowStorage.cpp
    feEigenProblem.cpp
    feExporter.cpp
    feGraphics.cpp
    feLinearSystem.cpp
    feLinearSystemMklPardiso.cpp
    feLinearSystemPETSc.cpp
    feMatrixInterface.cpp
    feMesh.cpp
    feMeshRead.cpp
    feMessage.cpp
    feMetric.cpp
    feMetricAnalyticP2.cpp
    feMetricComputation.cpp
    feMetricInterpolation.cpp
    feMetricGradation.cpp
    feMetricTools.cpp
    feMinimizeInterpolationError.cpp
    feNG.cpp
    feNonLinearSolver.cpp
    feNorm.cpp
    feNumber.cpp
    feNumeric.cpp
    feOptionsParser.cpp
    fePolynomial.cpp
    feRecovery.cpp
    feNewRecovery.cpp
    feQuadrature.cpp
    feQuadratureTet.cpp
    feQuadratureTri.cpp
    feSolution.cpp
    feSolutionContainer.cpp
    feSpace.cpp
    feSpace_1D.cpp
    feSpace_2D.cpp
    feSpace_3D.cpp
    feSysElm.cpp
    feSysElm_1D.cpp
    feSysElm_2D.cpp
    feTimeIntegration.cpp
    feTriangle.cpp
    feTetrahedron.cpp
    feVectorSpace.cpp
    feVectorSysElm.cpp
    GenericSolver.cpp
    geodesic.cpp
    glfem.cpp
    TransientAdapter.cpp)

#Add Gmsh files for STensor3
list(APPEND feNG_SRC ${CMAKE_SOURCE_DIR}/contrib/gmsh/fullMatrix.cpp
    # ${CMAKE_SOURCE_DIR}/contrib/gmsh/poly34.cpp
    # ${CMAKE_SOURCE_DIR}/contrib/gmsh/Numeric.cpp
    ${CMAKE_SOURCE_DIR}/contrib/gmsh/STensor3.cpp)

set(feNG_API
    CHNS_Solver.h
    DiffusionSolver.h
    ellipseToolbox.h
    feAPI.h
    feBilinearForm.h
    feBoundaryConditions.h
    feCncGeo.h
    feColoring.h
    feCompressedRowStorage.h
    feConstants.h
    feEdge.h
    feEigenProblem.h
    feExporter.h
    feGraphics.h
    feFunction.h
    feLinearSystem.h
    feMatrixInterface.h
    feMesh.h
    feMessage.h
    feMetric.h
    feMetricTools.h
    feMinimizeInterpolationError.h
    feNG.h
    feNonLinearSolver.cpp
    feNorm.h
    feNumber.h
    feNumeric.h
    feOptionsParser.h
    fePolynomial.h
    feQuadrature.h
    feRecovery.h
    feNewRecovery.h
    feSolution.h
    feSolutionContainer.h
    feSpace.h
    feSysElm.h
    feTimeIntegration.h
    feTriangle.h
    feTetrahedron.h
    feVertex.h
    GenericSolver.h
    geodesic.h
    glfem.h
    TransientAdapter.h)

addLibraryMacro(feNG "${feNG_SRC}" "${feNG_API}")
