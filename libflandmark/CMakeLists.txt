find_package( OpenCV REQUIRED )
include_directories(${OPENCV_INCLUDE_DIR})

if(${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 cxcore cv cvaux)
	elseif (WIN32)
		set(CV_LIBS_1 cxxore${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} cv${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} cvaux${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
else(${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 opencv_core opencv_imgproc)
	elseif (WIN32)
		set(CV_LIBS_1 opencv_core${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} opencv_imgproc${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
endif(${OpenCV_VERSION_MINOR} LESS 3)

add_library(flandmark_static STATIC flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
target_link_libraries(flandmark_static ${CV_LIBS_1})
if(CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(flandmark_static PROPERTIES COMPILE_FLAGS -fPIC)
endif(CMAKE_COMPILER_IS_GNUCC)
set_property(TARGET flandmark_static PROPERTY COMPILE_DEFINITIONS FLANDMARK_STATIC)

add_library(flandmark_shared SHARED flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
target_link_libraries(flandmark_shared ${CV_LIBS_1})

#setup Config.cmake
SET(FLANDMARK_BASE_DIR "${PROJECT_SOURCE_DIR}/libflandmark")
set(FLANDMARK_BINARY_DIR "${PROJECT_BINARY_DIR}/libflandmark")
configure_file(flandmarkConfig.cmake.in
  "${PROJECT_BINARY_DIR}/libflandmark/flandmarkConfig.cmake" @ONLY)

export(PACKAGE flandmark)
