# 
# These tests check the functionalities of the core library and 
# do not require a linear solver
# 
if(BUILD_TESTING)

	include(GoogleTest)

	# Add all .cpp files in this directory as test files
	file( GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/tests/withoutLinearSolver/*.cpp )
	foreach( file ${EXECUTABLES})
		get_filename_component(barename ${file} NAME)
	    string( REPLACE ".cpp" "" name ${barename})
	    message(STATUS "Adding test file ${name}")
	    feNG_add_executable( ${name} ${barename} GTest::gtest_main)
	    gtest_discover_tests(${name})
	endforeach( file ${EXECUTABLES} )

endif()