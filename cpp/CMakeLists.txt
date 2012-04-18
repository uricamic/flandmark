# create static library libflandmark_static
INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/cpp/libflandmark")
ADD_SUBDIRECTORY (libflandmark) 

# create mex-files flandmark_load_model and flandmark_detecor
IF (BUILD_MATLAB_BINDINGS)
  INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/cpp/mex")
  ADD_SUBDIRECTORY (mex)
ENDIF()

# simple example
FIND_PACKAGE( OpenCV REQUIRED )
INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIR})

set(${PROJECT_NAME}_simple_example_srcs simple_example.cpp)
add_executable(${PROJECT_NAME}_simple_example ${${PROJECT_NAME}_simple_example_srcs})
target_link_libraries(${PROJECT_NAME}_simple_example flandmark_static opencv_highgui)

set(${PROJECT_NAME}_example1_srcs example1.cpp)
add_executable(${PROJECT_NAME}_1 ${${PROJECT_NAME}_example1_srcs})
target_link_libraries(${PROJECT_NAME}_1 flandmark_static opencv_highgui opencv_core opencv_objdetect)

set(${PROJECT_NAME}_example2_srcs example2.cpp)
add_executable(${PROJECT_NAME}_2 ${${PROJECT_NAME}_example2_srcs})
target_link_libraries(${PROJECT_NAME}_2 flandmark_static opencv_highgui opencv_core opencv_objdetect)
