include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

if(ENABLE_STATIONARY)
	add_subdirectory(S)
endif()

if(ENABLE_TRANSIENT)
	add_subdirectory(T)
endif()

if(ENABLE_GRAPHICS)
	add_subdirectory(G)
endif()

if(ENABLE_INTERPOLATION)
	add_subdirectory(I)
endif()