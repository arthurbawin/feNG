cmake_minimum_required(VERSION 3.3)
project(feNG)

###################################
# Options
###################################
option(ENABLE_OPENMP        "Enable OpenMP"               OFF)
option(ENABLE_MPI           "Enable MPI"                  ON )
option(ENABLE_BOOST         "Enable Boost"                OFF)
option(ENABLE_SUNDIALS      "Enable SUNDIALS"             OFF)
option(ENABLE_SOPLEX        "Enable SOPLEX"               OFF)
option(ENABLE_PETSC         "Enable PETSc"                ON )
option(ENABLE_MKL           "Enable MKL"                  OFF)
option(ENABLE_SANITY_CHECKS "Enable sanity checks"        OFF )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

#### Organize outputs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

###################################
####  Compiler-dependent flags ####
###################################

message(STATUS "========== INFO : Compilateur C   = ${CMAKE_C_COMPILER}")
message(STATUS "========== INFO : Compilateur C++ = ${CMAKE_CXX_COMPILER}")

#######################
#  GNU compiler flags #
#######################
if ((NOT MSVC) AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")

  #### Standard C++11
  enable_language(CXX)
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)

  #### Options de compilation
  add_compile_options("-ansi;-m64")

  # Crash à l'exécution quand on compile en O3 ?
  #add_compile_options("-O3;-ansi;-m64")
  add_compile_options("$<$<CONFIG:RELEASE>:-ansi;-m64;-O3;-fast>"
                      "$<$<CONFIG:RELWITHDEBINFO>:-Wall;-ansi;-m64;-O3;-fast;-g>"
                      "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-ansi;-m64;-O0;-g>"
                      #"-Wall;-Wextra;-O3;-ansi;-m64"
			)

  #### Sanity checks
  if(ENABLE_SANITY_CHECKS)
    add_compile_options("-Wall"
                        "-Wextra"
                        "-Wno-sign-compare"
                        "-Wno-unused-parameter"
                        "-fno-omit-frame-pointer"
                        "-Wno-missing-field-initializers"
                        "-Wpointer-arith"
                        "-Wcast-align"
                        # "-Wcast-qual"
                        "-Wstrict-aliasing"
                        "-Wpointer-arith"
                        "-Winit-self"
                        # "-Wredundant-decls"
                        # "-fsanitize=address"
                        # "-fsanitize=leak"
                        # "-fsanitize=undefined"
                        )
    
    # set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
    #   -fsanitize=address -fsanitize=leak -fsanitize=undefined")
  endif()
  
########################
# Intel compiler flags #
########################
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")

  #### Options de compilation
  add_compile_options("$<$<CONFIG:RELEASE>:-ansi;-m64;-O3;-fast>"
                      "$<$<CONFIG:RELWITHDEBINFO>:-Wall;-ansi;-m64;-O3;-fast;-g>"
                      "$<$<CONFIG:DEBUG>:-Wall;-Wextra;-ansi;-m64;-O0;-g;-debug full>"
		      "-Wall;-Wextra;-O3;-ansi;-m64")

	      add_compile_options("-no-inline-min-size"
	              "-no-inline-max-size"
	              "-no-inline-max-total-size")

  #### Standard C++11
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

  if(COMPILER_SUPPORTS_CXX11)
    message(STATUS "Ajout du flag -std=c++11\n")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    add_compile_options("-std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    message(STATUS "Ajout du flag -std=c++0x\n")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    add_compile_options("-std=c++0x")
  else()
      message(WARNING "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()

#######################
# MSVC compiler flags #
#######################
elseif (MSVC)

  message(STATUS "==================================================")
  message(STATUS "Options de compilation MSVC encore experimentales ")
  message(STATUS "==================================================")

  #### Options de compilation
  add_compile_options("$<$<OR:$<CONFIG:RELEASE>,
                              $<CONFIG:RELWITHDEBINFO>>:-O2;-fp:strict;-GL>"
                      "-MP")

##########################
# Other compiler warning #
##########################
else()
  message(WARNING "Unsupported compiler")
endif()

############################################
#### End of Platform dependent settings ####
############################################

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_definitions(-DDEBUG)
endif()

#### Initialize list of external libraries
set(feNG_EXTERNAL_LIBS "")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

############################################
#  OpenMP
############################################
if(ENABLE_OPENMP)
  if(MSVC AND CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(OMP_LIB_PATH "C:/Program\ Files/LLVM/lib")
    add_library(omp_lib STATIC IMPORTED) # or STATIC instead of SHARED
    set_target_properties(omp_lib PROPERTIES
      IMPORTED_LOCATION "${OMP_LIB_PATH}/libomp.lib"
    )
    list(APPEND HXT_EXTERNAL_LIBS omp_lib)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    add_definitions("-DHAVE_OPENMP")
  else()
    FIND_PACKAGE(OpenMP 4)
    if(OPENMP_FOUND)
    message(STATUS "OpenMP trouvé\n")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
      add_definitions("-DHAVE_OPENMP")
    else()
      message(STATUS "OpenMP non trouvé\n")
    endif()
  endif()
endif()

############################################
#  MPI
############################################
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

############################################
#  MKL
############################################
if(ENABLE_MKL)
  #include(CheckLibraryExists)
  #include(LibraryUtils)
  #include(CheckLinkerFlag)
  if(ENABLE_OPENMP)
    set(MKL_USE_parallel 1)
  endif()
  find_package(MKL REQUIRED)
  if(MKL_FOUND)
    message(STATUS "MKL trouve ! (-:")
    include_directories(${MKL_INCLUDE_DIRS})
    list(APPEND EF_EXTERNAL_LIBS ${EF_EXTERNAL_LIBS} ${MKL_LIBRARIES})
    add_definitions("-DHAVE_MKL")
    else()
    message(STATUS "MKL non trouve )-:")
  endif()
endif()

############################################
#  PETSc
############################################
if(ENABLE_PETSC)
  if(NOT ENABLE_MPI)
    MESSAGE(FATAL_ERROR "PETSc requires MPI : please turn on ENABLE_MPI.")
  endif()
  MESSAGE(STATUS "petsc_found : ${PETSC_FOUND}")
  unset(PETSC_FOUND CACHE)
  MESSAGE(STATUS "petsc_found : ${PETSC_FOUND}")
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig")
  MESSAGE(STATUS "petsc dir : ${PETSC_DIR}")
  MESSAGE(STATUS "petsc arch : ${PETSC_ARCH}")
  MESSAGE(STATUS "pkg config path : $ENV{PKG_CONFIG_PATH}")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(PETSC REQUIRED petsc)
  include_directories(${PETSC_INCLUDE_DIRS})
  MESSAGE(STATUS "petsc infos : ${PETSC_FOUND}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_LIBRARIES}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_LINK_LIBRARIES}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_LIBRARY_DIRS}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_LDFLAGS}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_LDFLAGS_OTHER}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_INCLUDE_DIRS}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_CFLAGS}")
  # MESSAGE(STATUS "petsc infos : ${PETSC_CFLAGS_OTHER}")
  list(APPEND feNG_EXTERNAL_LIBS ${PETSC_LINK_LIBRARIES} ${PETSC_LDFLAGS})
  add_definitions("-DHAVE_PETSC")
  MESSAGE(STATUS "PETSc found (-:")
endif()

############################################
#  BOOST
############################################
if(ENABLE_BOOST)
  find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${Boost_LIBRARIES})
  add_definitions("-DHAVE_BOOST")
  MESSAGE(STATUS "BOOST found (-:")
endif()

############################################
#  SUNDIALS
############################################
if(ENABLE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
  include_directories(${SUNDIALS_INCLUDE_DIRS})
  list(APPEND feNG_EXTERNAL_LIBS ${SUNDIALS_LIBRARIES})
  add_definitions("-DHAVE_SUNDIALS")
  MESSAGE(STATUS "SUNDIALS found (-:")
endif()

############################################
#  SoPlex (algorithme du simplexe)
############################################
if(ENABLE_SOPLEX)
  list(APPEND feNG_EXTERNAL_LIBS -lsoplex)
  add_definitions("-DHAVE_SOPLEX")
  MESSAGE(STATUS "SOPLEX found (-:")
endif()

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

############################################
#  Macro pour ajouter un executable
############################################
macro(feNG_add_executable EXENAME SRC LIBS)
  add_executable(${EXENAME} ${SRC})
  target_link_libraries(${EXENAME} ${LIBS} feNG ${feNG_EXTERNAL_LIBS})
  install(TARGETS ${EXENAME} DESTINATION bin)
endmacro()

MESSAGE(STATUS "External libs for feNG : ${feNG_EXTERNAL_LIBS}")

#### Sous-repertoires
add_subdirectory(src)
add_subdirectory(exe)