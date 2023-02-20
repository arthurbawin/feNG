include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

file( GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/exe/*.cpp)
foreach( file ${EXECUTABLES})
	get_filename_component(barename ${file} NAME)
    string( REPLACE ".cpp" "" name ${barename})
    message(STATUS "Adding executable ${name}")
    feNG_add_executable( ${name} ${barename} "")
endforeach( file ${EXECUTABLES} )