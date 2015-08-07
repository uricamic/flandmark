# Compile examples

include_directories("${PROJECT_SOURCE_DIR}/libflandmark")

SET( OpenCV "E:/Program Files/opencv")
SET( OPENCV_INCLUDE_DIR "E:/Program Files/opencv/build/include")
include_directories(${OPENCV_INCLUDE_DIR})

set(CV_LIBS_1 opencv_highgui2410)
set(CV_LIBS_2 opencv_core2410 opencv_objdetect2410)
	
set(${PROJECT_NAME}_simple_example_srcs simple_example.cpp)
add_executable(${PROJECT_NAME}_simple_example ${${PROJECT_NAME}_simple_example_srcs})
target_link_libraries(${PROJECT_NAME}_simple_example flandmark_static ${CV_LIBS_1})

set(${PROJECT_NAME}_example1_srcs example1.cpp)
add_executable(${PROJECT_NAME}_1 ${${PROJECT_NAME}_example1_srcs})
target_link_libraries(${PROJECT_NAME}_1 flandmark_static ${CV_LIBS_1} ${CV_LIBS_2})

set(${PROJECT_NAME}_example2_srcs example2.cpp)
add_executable(${PROJECT_NAME}_2 ${${PROJECT_NAME}_example2_srcs})
target_link_libraries(${PROJECT_NAME}_2 flandmark_static ${CV_LIBS_1} ${CV_LIBS_2})

