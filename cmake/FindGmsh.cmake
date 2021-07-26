##############################################################
# Try to find Gmsh                                           #
#                                                            #
# Once done this will define:                                #
#  GMSH_FOUND   - system has Gmsh                            #
#  GMSH_INC     - Gmsh include directory (static or dynamic) #
#  GMSH_LIB     - Gmsh library                               #
#                                                            #
# Usage:                                                     #
#  find_package(Gmsh)                                        #
#                                                            #
# Setting these changes the behavior of the search           #
#  GMSH_INC     - Gmsh include directory                     #
#  GMSH_LIB     - Gmsh library path (static or dynamic)      #
##############################################################

## Try to set GMSH_LIB and GMSH_INC from environment variables ##
#################################################################
if(NOT DEFINED GMSH_LIB)
  set(GMSH_LIB $ENV{GMSH_LIB})
endif()
if(NOT DEFINED GMSH_INC)
  set(GMSH_INC $ENV{GMSH_INC})
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
  "Gmsh could not be found: be sure to set GMSH_LIB and GMSH_INC in your environment variables"
  GMSH_LIB GMSH_INC)
