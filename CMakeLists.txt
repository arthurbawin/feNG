cmake_minimum_required(VERSION 3.3)
project(feNG)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_METRIC        "Enable metric tensor computations" ON )
option(ENABLE_OPENMP        "Enable OpenMP"                     ON )
option(ENABLE_MPI           "Enable MPI"                        ON )
option(ENABLE_PETSC         "Enable PETSc"                      ON )
option(ENABLE_SLEPC         "Enable SLEPc"                      OFF)
option(ENABLE_GMSH          "Enable Gmsh"                       OFF)
option(ENABLE_MKL           "Enable MKL"                        OFF)
option(ENABLE_METIS         "Enable Metis"                      OFF)
option(ENABLE_SOPLEX        "Enable SoPlex"                     OFF)
option(ENABLE_PROFILER      "Compile with -lprofiler"           OFF)
option(ENABLE_GLFW          "Enable basic GLFW graphics"        OFF)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

  # #### Options de compilation
  # add_compile_options("-ansi;-m64")

  # add_compile_options("$<$<CONFIG:RELEASE>:-ansi;-m64;-O3;-fast>"
  #                     "$<$<CONFIG:RELWITHDEBINFO>:-Wall;-ansi;-m64;-O3;-fast;-g>"
  #                     "$<$<CONFIG:DEBUG>:-ansi;-m64;-O0;-g;-ggdb3>"
  #                     # "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-ansi;-m64;-O0;-g;-ggdb3>"
  #                     #"-Wall;-Wextra;-O3;-ansi;-m64"
  #                     )

# add_compile_options(-pedantic)
# add_compile_options(-pedantic-errors)
# add_compile_options(-Wextra)
add_compile_options(-Wreturn-type)
add_compile_options(-Wunused-variable)
add_compile_options(-Wdisabled-optimization)
# add_compile_options(-Werror)

############################################
#  Macro pour definir une bibliotheque
############################################
macro(feNG_add_library LIBNAME SRC API_H)
  foreach(HFILE ${API_H})
    list(APPEND feNG_API_H "${CMAKE_CURRENT_SOURCE_DIR}/${HFILE}")
  endforeach()
  SET(feNG_API_H ${feNG_API_H} PARENT_SCOPE)
  add_library(${LIBNAME} ${SRC})
  install(TARGETS ${LIBNAME} DESTINATION lib)
endmacro()


#  Macro pour ajouter un executable
macro(feNG_add_executable EXENAME SRC LIBS)
  add_executable(${EXENAME} ${SRC})
  target_link_libraries(${EXENAME} ${LIBS} feNG ${feNG_EXTERNAL_LIBS})
  install(TARGETS ${EXENAME} DESTINATION bin)
endmacro()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DHAVE_DEBUG)
endif()

# Initialize list of external libraries
set(feNG_EXTERNAL_LIBS "")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# It seems that CMake only defines -DNDEBUG when it is NOT in debug mode
# Add also a -DFENG_DEBUG definition to avoid double negative conditions
# -DDEBUG conflicts with definitions from other libraries such as SoPlex
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DFENG_DEBUG)
endif()

#  OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    MESSAGE(STATUS "OpenMP found (-:")
    list(APPEND feNG_EXTERNAL_LIBS OpenMP::OpenMP_CXX)
    add_definitions("-DHAVE_OMP")
    add_compile_options(-fopenmp)
  endif()
endif()

#  MPI
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  MESSAGE(STATUS "MPI infos: ${MPI_CXX_INCLUDE_PATH}")
  MESSAGE(STATUS "MPI infos: ${MPI_LIBRARIES}")
  MESSAGE(STATUS "MPI infos: ${MPI_CXX_LIBRARIES}")
  MESSAGE(STATUS "MPI infos: ${MPI_CXX_LINK_FLAGS}")
  MESSAGE(STATUS "MPI infos: ${MPI_CFLAGS_OTHER}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  list(APPEND feNG_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  add_definitions("-DHAVE_MPI")
  MESSAGE(STATUS "MPI found (-:")
endif()

#  MKL
if(ENABLE_MKL)
  # if(NOT ENABLE_MPI)
  #   MESSAGE(FATAL_ERROR "MKL requires MPI : please turn on ENABLE_MPI.")
  # endif()
  # include(CheckLibraryExists)
  # include(LibraryUtils)
  # include(CheckLinkerFlag)
  if(ENABLE_OPENMP)
    set(MKL_USE_parallel 1)
    set(MKL_OPENMP "gomp")
  endif()
  find_package(MKL REQUIRED)
  message(STATUS "Found MKL_FOUND:\n${MKL_FOUND} ")
  message(STATUS "Found MKL_ROOT_DIR:\n${MKL_ROOT_DIR} ")
  message(STATUS "Found MKL_INCLUDE_DIR:\n${MKL_INCLUDE_DIR} ")
  message(STATUS "Found MKL_LIBRARIES:\n${MKL_LIBRARIES} ")
  message(STATUS "Found MKL_LIBRARY_DIR:\n${MKL_LIBRARY_DIR} ")
  include_directories(${MKL_INCLUDE_DIR})
  list(APPEND feNG_EXTERNAL_LIBS ${feNG_EXTERNAL_LIBS} ${MKL_LIBRARIES})
  add_definitions("-DHAVE_MKL")
  message(STATUS "MKL found ! (-:")
endif()

#  PETSc
if(ENABLE_PETSC)
  if(NOT ENABLE_MPI)
    MESSAGE(FATAL_ERROR "PETSc requires MPI : please turn on ENABLE_MPI.")
  endif()
  string(REGEX REPLACE "/+$" "" PETSC_DIR "${PETSC_DIR}")
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig")
  # MESSAGE(STATUS "PETSC_DIR            : ${PETSC_DIR}")
  # MESSAGE(STATUS "PETSC_ARCH           : ${PETSC_ARCH}")
  # MESSAGE(STATUS "PKG_CONFIG_PATH      : $ENV{PKG_CONFIG_PATH}")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(PETSC REQUIRED PETSc)
  include_directories(${PETSC_INCLUDE_DIRS})
  # MESSAGE(STATUS "PETSC_FOUND          : ${PETSC_FOUND}")
  # MESSAGE(STATUS "PETSC_LIBRARIES      : ${PETSC_LIBRARIES}")
  # MESSAGE(STATUS "PETSC_LINK_LIBRARIES : ${PETSC_LINK_LIBRARIES}")
  # MESSAGE(STATUS "PETSC_LIBRARY_DIRS   : ${PETSC_LIBRARY_DIRS}")
  # MESSAGE(STATUS "PETSC_LDFLAGS        : ${PETSC_LDFLAGS}")
  # MESSAGE(STATUS "PETSC_LDFLAGS_OTHER  : ${PETSC_LDFLAGS_OTHER}")
  # MESSAGE(STATUS "PETSC_INCLUDE_DIRS   : ${PETSC_INCLUDE_DIRS}")
  # MESSAGE(STATUS "PETSC_CFLAGS         : ${PETSC_CFLAGS}")
  # MESSAGE(STATUS "PETSC_CFLAGS_OTHER   : ${PETSC_CFLAGS_OTHER}")
  list(APPEND feNG_EXTERNAL_LIBS ${PETSC_LINK_LIBRARIES} ${PETSC_LDFLAGS})
  add_definitions("-DHAVE_PETSC")
  MESSAGE(STATUS "PETSc found (-:")
endif()

#  SLEPc (Eigenvalue solver based on PETSc)
if(ENABLE_SLEPC)
  if(NOT ENABLE_PETSC)
    MESSAGE(FATAL_ERROR "SLEPc requires PETSc : please turn on ENABLE_PETSC.")
  endif()
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${SLEPC_DIR}/${SLEPC_ARCH}/lib/pkgconfig")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SLEPC REQUIRED slepc)
  include_directories(${SLEPC_INCLUDE_DIRS})
  # MESSAGE(STATUS "PETSC_FOUND          : ${PETSC_FOUND}")
  # MESSAGE(STATUS "PETSC_LIBRARIES      : ${PETSC_LIBRARIES}")
  # MESSAGE(STATUS "PETSC_LINK_LIBRARIES : ${PETSC_LINK_LIBRARIES}")
  # MESSAGE(STATUS "PETSC_LIBRARY_DIRS   : ${PETSC_LIBRARY_DIRS}")
  # MESSAGE(STATUS "PETSC_LDFLAGS        : ${PETSC_LDFLAGS}")
  # MESSAGE(STATUS "PETSC_LDFLAGS_OTHER  : ${PETSC_LDFLAGS_OTHER}")
  # MESSAGE(STATUS "PETSC_INCLUDE_DIRS   : ${PETSC_INCLUDE_DIRS}")
  # MESSAGE(STATUS "PETSC_CFLAGS         : ${PETSC_CFLAGS}")
  # MESSAGE(STATUS "PETSC_CFLAGS_OTHER   : ${PETSC_CFLAGS_OTHER}")
  list(APPEND feNG_EXTERNAL_LIBS ${SLEPC_LINK_LIBRARIES} ${SLEPC_LDFLAGS})
  add_definitions("-DHAVE_SLEPC")
  MESSAGE(STATUS "SLEPc found (-:")
endif()

#  Metis
if(ENABLE_METIS)
  find_package(METIS)
  MESSAGE(STATUS "METIS_FOUND          : ${METIS_FOUND}")
  MESSAGE(STATUS "METIS_INCLUDE_DIRS   : ${METIS_INCLUDE_DIRS}")
  MESSAGE(STATUS "METIS_LIBRARY_DIRS   : ${METIS_LIBRARY_DIRS}")
  MESSAGE(STATUS "METIS_LIBRARIES      : ${METIS_LIBRARIES}")
  include_directories(${METIS_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${METIS_LIBRARIES})
  add_definitions("-DHAVE_METIS")
  MESSAGE(STATUS "METIS found (-:")
endif()

#  SoPlex (optimization library for aniso and curved mesh adaptation)
if(ENABLE_SOPLEX)
  find_package(SOPLEX)
  MESSAGE(STATUS "SOPLEX_FOUND          : ${SOPLEX_FOUND}")
  MESSAGE(STATUS "SOPLEX_INCLUDE_DIRS   : ${SOPLEX_INCLUDE_DIRS}")
  MESSAGE(STATUS "SOPLEX_PIC_LIBRARIES  : ${SOPLEX_PIC_LIBRARIES}")
  MESSAGE(STATUS "SOPLEX_LIBRARIES      : ${SOPLEX_LIBRARIES}")
  include_directories(${SOPLEX_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${SOPLEX_LIBRARIES} ${SOPLEX_PIC_LIBRARIES})
  # GMP is required by SoPlex to solve with higher precision
  # Tried to compile SoPlex with GMP=false but it still seems required
  add_compile_options("-lgmp")
  add_definitions("-DHAVE_SOPLEX")
  MESSAGE(STATUS "SoPlex found (-:")
endif()

#  Eigen
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/Eigen)
add_definitions("-DHAVE_EIGEN")

#  Root-finder (https://github.com/ZhepeiWang/Root-Finder)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/root-finder)

#  Gmsh
if(ENABLE_GMSH)
  find_package(Gmsh)
  set(GMSH_LIB "${GMSH_LIB}" CACHE FILEPATH "Gmsh library")
  set(GMSH_INC "${GMSH_INC}" CACHE PATH     "Gmsh headers")
  include_directories(${GMSH_INC})
  list(APPEND feNG_EXTERNAL_LIBS ${GMSH_LIB})
  add_definitions("-DHAVE_GMSH")
  MESSAGE(STATUS "Gmsh found (-:")
endif()

# Gmsh classes for metric tensors and rtree 
# See also src/CmakeLists.txt for the .cpp files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/gmsh)

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

MESSAGE(STATUS "External libs for feNG : ${feNG_EXTERNAL_LIBS}")

add_subdirectory(src)
add_subdirectory(exe)
