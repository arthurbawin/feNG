cmake_minimum_required(VERSION 3.3)
project(feNG)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_OPENMP        "Enable OpenMP"               OFF)
option(ENABLE_MPI           "Enable MPI"                  ON )
option(ENABLE_PETSC         "Enable PETSc"                ON )
option(ENABLE_MKL           "Enable MKL"                  OFF)

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

add_compile_options(-pedantic)
add_compile_options(-pedantic-errors)
add_compile_options(-Wextra)
add_compile_options(-Wreturn-type)
add_compile_options(-Wunused-variable)
add_compile_options(-Wdisabled-optimization)
add_compile_options(-fext-numeric-literals)

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
  add_definitions(-DDEBUG)
endif()

# Initialize list of external libraries
set(feNG_EXTERNAL_LIBS "")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#  MPI
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  # MESSAGE(STATUS "MPI infos: ${MPI_CXX_INCLUDE_PATH}")
  # MESSAGE(STATUS "MPI infos: ${MPI_LIBRARIES}")
  # MESSAGE(STATUS "MPI infos: ${MPI_CXX_LIBRARIES}")
  # MESSAGE(STATUS "MPI infos: ${MPI_CXX_LINK_FLAGS}")
  # MESSAGE(STATUS "MPI infos: ${MPI_CFLAGS_OTHER}")
  include_directories(${MPI_CXX_INCLUDE_PATH})
  list(APPEND feNG_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
  add_definitions("-DHAVE_MPI")
  MESSAGE(STATUS "MPI found (-:")
endif()

#  MKL
if(ENABLE_MKL)
  #include(CheckLibraryExists)
  #include(LibraryUtils)
  #include(CheckLinkerFlag)
  # if(ENABLE_OPENMP)
  #   set(MKL_USE_parallel 1)
  # endif()
  find_package(MKL REQUIRED)
  include_directories(${MKL_INCLUDE_DIRS})
  list(APPEND EF_EXTERNAL_LIBS ${EF_EXTERNAL_LIBS} ${MKL_LIBRARIES})
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

#  Eigen
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/Eigen)

#  Gmsh classes for metric tensors and rtree
# See also src/CmakeLists.txt for the .cpp files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/contrib/gmsh)

MESSAGE(STATUS "External libs for feNG : ${feNG_EXTERNAL_LIBS}")

add_subdirectory(src)
add_subdirectory(exe)