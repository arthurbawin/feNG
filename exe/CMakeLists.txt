include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add all .cpp files in the directory
file( GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/exe/*.cpp)
foreach(file ${EXECUTABLES})
	get_filename_component(barename ${file} NAME)
    string(REPLACE ".cpp" "" name ${barename})
    message(STATUS "Adding executable ${name}")
    addExecutableMacro( ${name} ${barename} "")
endforeach( file ${EXECUTABLES} )