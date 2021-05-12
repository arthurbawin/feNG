include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

feNG_add_executable(testStationnaire testStationnaire.cpp "")
feNG_add_executable(testBDF2 testBDF2.cpp "")
feNG_add_executable(testDC3 testDC3.cpp "")
feNG_add_executable(testCSR testCSR.cpp "")

#target_include_directories(feNGt PUBLIC ${MKL_INCLUDE_DIRS})
#target_link_libraries(feNGs ${MKL_LIBRARIES})