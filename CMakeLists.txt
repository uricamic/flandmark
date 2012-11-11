cmake_minimum_required(VERSION 2.6)

project(flandmark)

# The version number.
set(flandmark_VERSION_MAJOR 1)
set(flandmark_VERSION_MINOR 7)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_MATLAB_BINDINGS "Build MATLAB bindings" ON)
option(BUILD_CPP_EXAMPLES "Build CPP examples" ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")

# build libflandmark shared and static libraries
add_subdirectory(libflandmark)
# build cpp examples
if (BUILD_CPP_EXAMPLES)
	add_subdirectory(examples)
endif(BUILD_CPP_EXAMPLES)
# build MATLAB mex-files
if(BUILD_MATLAB_BINDINGS)
	add_subdirectory(matlab_toolbox/mex)
endif(BUILD_MATLAB_BINDINGS)

# copy some important files to the binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/flandmark_model.dat ${CMAKE_CURRENT_BINARY_DIR}/examples/flandmark_model.dat COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/face.jpg ${CMAKE_CURRENT_BINARY_DIR}/examples/face.jpg COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/groupphoto.jpg ${CMAKE_CURRENT_BINARY_DIR}/examples/groupphoto.jpg COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/seq_bruges04_300frames.avi ${CMAKE_CURRENT_BINARY_DIR}/examples/seq_bruges04_300frames.avi COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/haarcascade_frontalface_alt.xml ${CMAKE_CURRENT_BINARY_DIR}/examples/haarcascade_frontalface_alt.xml COPYONLY)

