# Setting this will only affect the folders down from the current one
set(CMAKE_CXX_STANDARD 17)
add_subdirectory(googletest)

if(BUILD_TESTING)

	include(GoogleTest)

	# Add all .cpp files in this directory as test files
	file( GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/tests/*.cpp )
	foreach( file ${EXECUTABLES})
		get_filename_component(barename ${file} NAME)
	    string( REPLACE ".cpp" "" name ${barename})
	    message(STATUS "Adding test file ${name}")
	    feNG_add_executable( ${name} ${barename} GTest::gtest_main)
	    gtest_discover_tests(${name})
	endforeach( file ${EXECUTABLES} )

	# add_executable(my_test my_test.cpp)
	# target_link_libraries(my_test PRIVATE feNG GTest::gtest_main)

	# include(GoogleTest)
	# gtest_discover_tests(my_test)
endif()