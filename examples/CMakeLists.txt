# Compile examples

include_directories("${PROJECT_SOURCE_DIR}/libflandmark")

find_package( OpenCV REQUIRED )
include_directories(${OPENCV_INCLUDE_DIR})

if (${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 highgui)
		set(CV_LIBS_2 cxcore cv cvaux)
	elseif (WIN32)
		set(CV_LIBS_1 highgui${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
		set(CV_LIBS_2 cxcore${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} cv${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} cvaux${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
else(${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 opencv_highgui)
		set(CV_LIBS_2 opencv_core opencv_objdetect)
	elseif (WIN32)
		set(CV_LIBS_1 opencv_highgui${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
		set(CV_LIBS_2 opencv_core${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} opencv_objdetect${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
endif(${OpenCV_VERSION_MINOR} LESS 3)

set(${PROJECT_NAME}_simple_example_srcs simple_example.cpp)
add_executable(${PROJECT_NAME}_simple_example ${${PROJECT_NAME}_simple_example_srcs})
target_link_libraries(${PROJECT_NAME}_simple_example flandmark_static ${CV_LIBS_1})

set(${PROJECT_NAME}_example1_srcs example1.cpp)
add_executable(${PROJECT_NAME}_1 ${${PROJECT_NAME}_example1_srcs})
target_link_libraries(${PROJECT_NAME}_1 flandmark_static ${CV_LIBS_1} ${CV_LIBS_2})

set(${PROJECT_NAME}_example2_srcs example2.cpp)
add_executable(${PROJECT_NAME}_2 ${${PROJECT_NAME}_example2_srcs})
target_link_libraries(${PROJECT_NAME}_2 flandmark_static ${CV_LIBS_1} ${CV_LIBS_2})

