cmake_minimum_required(VERSION 3.3)
project(feNG)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TESTING       "Enable testing"                    OFF)
option(ENABLE_OPENMP        "Enable OpenMP"                     OFF)
option(ENABLE_MPI           "Enable MPI"                        OFF)
option(ENABLE_PETSC         "Enable PETSc"                      OFF)
option(ENABLE_SLEPC         "Enable SLEPc"                      OFF)
option(ENABLE_GMSH          "Enable Gmsh"                       OFF)
option(ENABLE_CURVED_MESH   "Enable curved meshing with Gmsh"   OFF)
option(ENABLE_MKL           "Enable MKL"                        OFF)
option(ENABLE_METIS         "Enable Metis"                      OFF)
option(ENABLE_SOPLEX        "Enable SoPlex"                     OFF)
option(ENABLE_CERES         "Enable Ceres-Solver"               OFF)
option(ENABLE_ALGLIB        "Enable AlgLib"                     OFF)
option(ENABLE_PROFILER      "Compile with -lprofiler"           OFF)
option(ENABLE_GLFW          "Enable basic GLFW graphics"        OFF)

# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

  # #### Options de compilation

  # add_compile_options("$<$<CONFIG:RELEASE>:-ansi;-m64;-O3;-fast>"
  #                     "$<$<CONFIG:RELWITHDEBINFO>:-Wall;-ansi;-m64;-O3;-fast;-g>"
  #                     "$<$<CONFIG:DEBUG>:-ansi;-m64;-O0;-g;-ggdb3>"
  #                     # "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-ansi;-m64;-O0;-g;-ggdb3>"
  #                     #"-Wall;-Wextra;-O3;-ansi;-m64"
  #                     )

add_compile_options(-pedantic)
add_compile_options(-pedantic-errors)
add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-Werror=return-type)
add_compile_options(-Wno-error=bool-compare)
add_compile_options(-Wshadow)
add_compile_options(-Wdisabled-optimization)

# Files from contrib/:
#
# Eigen and "unsupported" Eigen routines
# It may be necessary to include Eigen as system headers (SYSTEM) to
# silence Eigen's warnings such as -Wshadow
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen)
add_definitions("-DHAVE_EIGEN")
#  Root-finder (https://github.com/ZhepeiWang/Root-Finder)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/root-finder)
# Gmsh classes for metric tensors and rtree 
# See also src/CmakeLists.txt for the .cpp files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/gmsh)

# Macro to define a library
# Adapted from the HXT project:
# https://gitlab.onelab.info/gmsh/gmsh/tree/master/contrib/hxt
macro(addLibraryMacro LIBNAME SRC API_H)
  foreach(HFILE ${API_H})
    list(APPEND feNG_API_H "${CMAKE_CURRENT_SOURCE_DIR}/${HFILE}")
  endforeach()
  SET(feNG_API_H ${feNG_API_H} PARENT_SCOPE)
  add_library(${LIBNAME} ${SRC})
  install(TARGETS ${LIBNAME} DESTINATION lib)
endmacro()

# Macro to define a library
macro(addExecutableMacro EXENAME SRC LIBS)
  add_executable(${EXENAME} ${SRC})
  target_link_libraries(${EXENAME} ${LIBS} feNG ${feNG_EXTERNAL_LIBS})
  install(TARGETS ${EXENAME} DESTINATION bin)
endmacro()

# Initialize list of external libraries
set(feNG_EXTERNAL_LIBS "")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# It seems that CMake only defines -DNDEBUG when it is NOT in debug mode
# Add also a -DFENG_DEBUG definition to avoid double negative conditions
# -DDEBUG conflicts with definitions from other libraries such as SoPlex
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DFENG_DEBUG)
endif()

# OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found (-:")
    list(APPEND feNG_EXTERNAL_LIBS OpenMP::OpenMP_CXX)
    add_compile_options(-fopenmp)
    add_definitions("-DHAVE_OMP")
  endif()
endif()

#  MPI
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "MPI_CXX_COMPILER    : ${MPI_CXX_COMPILER}    ")
  message(STATUS "MPI_CXX_INCLUDE_DIRS: ${MPI_CXX_INCLUDE_DIRS}")
  message(STATUS "MPI_CXX_LIBRARIES   : ${MPI_CXX_LIBRARIES}   ")
  message(STATUS "MPI_CXX_LINK_FLAGS  : ${MPI_CXX_LINK_FLAGS}  ")
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  add_definitions("-DHAVE_MPI")
  message(STATUS "MPI found (-:")
endif()

# MKL
# Starting with Intel oneAPI 2021.3, a MKLConfig.cmake file is included.
# This is the recommended way of linking with MKL, as it seems that no
# official FindMKL.cmake was ever provided by Intel.
# For older versions of MKL, we use the custom FindMKL.cmake included in this project,
# which unfortunately has different option names than MKLConfig.cmake.
if(ENABLE_MKL)

  set(MKL_INTERFACE "ilp64") # 64bits interface
  if(ENABLE_OPENMP)
    set(MKL_THREADING "gnu_thread") # Only with GNU OpenMP for now
    # FindMKL.cmake options
    set(MKL_USE_parallel 1)
    set(MKL_OPENMP "gomp")
  else()
    set(MKL_THREADING "sequential")
  endif()
  # Only with OpenMPI for now. Not sure how to detect the MPI implementation.
  set(MKL_MPI "openmpi")

  # Try to locate MKLConfig.cmake first
  find_package(MKL CONFIG)
  if(${MKL_FOUND})
    message(STATUS "MKL_IMPORTED_TARGETS: ${MKL_IMPORTED_TARGETS}")
    message(STATUS "MKL_INCLUDE: ${MKL_INCLUDE} ")
    add_compile_options($<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  else()
    # Fallback to FindMKL.cmake
    message(STATUS "Could not find MKLConfig.cmake, using cmake/FindMKL.cmake instead")
    find_package(MKL REQUIRED)
    message(STATUS "MKL_ROOT_DIR   : ${MKL_ROOT_DIR}    ")
    message(STATUS "MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR} ")
    message(STATUS "MKL_LIBRARIES  : ${MKL_LIBRARIES}   ")
    message(STATUS "MKL_LIBRARY_DIR: ${MKL_LIBRARY_DIR} ")
    include_directories(${MKL_INCLUDE_DIR})
    list(APPEND feNG_EXTERNAL_LIBS ${feNG_EXTERNAL_LIBS} ${MKL_LIBRARIES})
    # Cluster Pardiso
    if(ENABLE_MPI)
      # Temporary fix to add the cluster options
      add_compile_options(-m64 -I"${MKLROOT}/include")
      list(APPEND feNG_EXTERNAL_LIBS "-m64  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl")
      # list(APPEND feNG_EXTERNAL_LIBS "-m64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl")
    endif()
  endif()

  add_definitions("-DHAVE_MKL")
  message(STATUS "MKL found (-:")
endif()

#  PETSc
if(ENABLE_PETSC)
  if(NOT ENABLE_MPI)
    message(FATAL_ERROR "PETSc requires MPI : please turn on ENABLE_MPI.")
  endif()
  message(STATUS "PETSC_DIR            : ${PETSC_DIR}")
  message(STATUS "PETSC_ARCH           : ${PETSC_ARCH}")
  string(REGEX REPLACE "/+$" "" PETSC_DIR "${PETSC_DIR}")
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig")
  message(STATUS "PKG_CONFIG_PATH      : $ENV{PKG_CONFIG_PATH}")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(PETSC REQUIRED PETSc)
  include_directories(${PETSC_INCLUDE_DIRS})
  message(STATUS "PETSC_FOUND          : ${PETSC_FOUND}")
  message(STATUS "PETSC_LIBRARIES      : ${PETSC_LIBRARIES}")
  message(STATUS "PETSC_LINK_LIBRARIES : ${PETSC_LINK_LIBRARIES}")
  message(STATUS "PETSC_LIBRARY_DIRS   : ${PETSC_LIBRARY_DIRS}")
  message(STATUS "PETSC_LDFLAGS        : ${PETSC_LDFLAGS}")
  message(STATUS "PETSC_LDFLAGS_OTHER  : ${PETSC_LDFLAGS_OTHER}")
  message(STATUS "PETSC_INCLUDE_DIRS   : ${PETSC_INCLUDE_DIRS}")
  message(STATUS "PETSC_CFLAGS         : ${PETSC_CFLAGS}")
  message(STATUS "PETSC_CFLAGS_OTHER   : ${PETSC_CFLAGS_OTHER}")
  list(APPEND feNG_EXTERNAL_LIBS ${PETSC_LINK_LIBRARIES} ${PETSC_LDFLAGS})
  add_definitions("-DHAVE_PETSC")
  message(STATUS "PETSc found (-:")
endif()

#  SLEPc (Eigenvalue solver based on PETSc)
if(ENABLE_SLEPC)
  if(NOT ENABLE_PETSC)
    message(FATAL_ERROR "SLEPc requires PETSc : please turn on ENABLE_PETSC.")
  endif()
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${SLEPC_DIR}/${SLEPC_ARCH}/lib/pkgconfig")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SLEPC REQUIRED slepc)
  include_directories(${SLEPC_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${SLEPC_LINK_LIBRARIES} ${SLEPC_LDFLAGS})
  add_definitions("-DHAVE_SLEPC")
  message(STATUS "SLEPc found (-:")
endif()

#  Metis
if(ENABLE_METIS)
  find_package(METIS)
  message(STATUS "METIS_FOUND          : ${METIS_FOUND}")
  message(STATUS "METIS_INCLUDE_DIRS   : ${METIS_INCLUDE_DIRS}")
  message(STATUS "METIS_LIBRARY_DIRS   : ${METIS_LIBRARY_DIRS}")
  message(STATUS "METIS_LIBRARIES      : ${METIS_LIBRARIES}")
  include_directories(${METIS_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${METIS_LIBRARIES})
  add_definitions("-DHAVE_METIS")
  message(STATUS "METIS found (-:")
endif()

#  SoPlex (optimization library for aniso and curved mesh adaptation)
if(ENABLE_SOPLEX)
  find_package(SOPLEX)
  message(STATUS "SOPLEX_FOUND          : ${SOPLEX_FOUND}")
  message(STATUS "SOPLEX_INCLUDE_DIRS   : ${SOPLEX_INCLUDE_DIRS}")
  message(STATUS "SOPLEX_PIC_LIBRARIES  : ${SOPLEX_PIC_LIBRARIES}")
  message(STATUS "SOPLEX_LIBRARIES      : ${SOPLEX_LIBRARIES}")
  include_directories(${SOPLEX_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${SOPLEX_LIBRARIES} ${SOPLEX_PIC_LIBRARIES})
  # GMP is required by SoPlex to solve with higher precision
  # Tried to compile SoPlex with GMP=false but it still seems required
  add_compile_options("-lgmp")
  add_definitions("-DHAVE_SOPLEX")
  message(STATUS "SoPlex found (-:")
endif()

# Ceres (optimization library)
if(ENABLE_CERES)
  find_package(Ceres REQUIRED)
  message(STATUS "CERES_FOUND     : ${CERES_FOUND}")
  message(STATUS "CERES_VERSION   : ${CERES_VERSION}")
  message(STATUS "CERES_LIBRARIES : ${CERES_LIBRARIES}")
  if(CERES_FOUND)
    add_definitions("-DHAVE_CERES")
    list(APPEND feNG_EXTERNAL_LIBS ${CERES_LIBRARIES})
  endif()
endif()

#  Gmsh
if(ENABLE_GMSH)
  find_package(Gmsh)
  set(GMSH_LIB "${GMSH_LIB}" CACHE FILEPATH "Gmsh library")
  set(GMSH_INC "${GMSH_INC}" CACHE PATH     "Gmsh headers")
  include_directories(${GMSH_INC})
  list(APPEND feNG_EXTERNAL_LIBS ${GMSH_LIB})
  add_definitions("-DHAVE_GMSH")
  message(STATUS "Gmsh found (-:")
endif()

if(ENABLE_CURVED_MESH)
  add_definitions("-DGMSH_WITH_CURVED_MESHING")
endif()

# ALGlib, not needed and should be removed
if(ENABLE_ALGLIB)
  # AlgLib for unconstrained optimization
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/alglib)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/alglib)
  list(APPEND feNG_EXTERNAL_LIBS alglib)
endif()

if(ENABLE_PROFILER)
  # For profiling with Google Perftools
  list(APPEND feNG_EXTERNAL_LIBS "-Wl,-no-as-needed,-lprofiler")
endif()

if(ENABLE_GLFW)
  find_package(OpenGL REQUIRED)
  include_directories(${OPENGL_INCLUDE_DIR})
  add_subdirectory(contrib/glfw)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/glfw/include)
  list(APPEND feNG_EXTERNAL_LIBS glfw ${OPENGL_gl_LIBRARY})
  add_definitions("-DHAVE_GLFW")
endif()
 
message(STATUS "External libs for feNG : ${feNG_EXTERNAL_LIBS}")

add_subdirectory(src)
add_subdirectory(exe)

if(ENABLE_TESTING)
  # Testing inspired by
  # https://github.com/cpp-for-yourself/lectures-and-homeworks/blob/main/lectures/googletest.md
  # Update the submodules here (googletest)
  include(cmake/UpdateSubmodules.cmake)
  include(CTest)
  add_subdirectory(tests)
endif()