include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add each executable individually
feNG_add_executable(example1 example1.cpp "")
feNG_add_executable(example2 example2.cpp "")
feNG_add_executable(testRecovery testRecovery.cpp "")

# Or add all .cpp files in the directory
# file( GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/exe/*.cpp)
# foreach( file ${EXECUTABLES})
# 	get_filename_component(barename ${file} NAME)
#     string( REPLACE ".cpp" "" name ${barename})
#     message(STATUS "Adding executable ${name}")
#     feNG_add_executable( ${name} ${barename} "")
# endforeach( file ${EXECUTABLES} )