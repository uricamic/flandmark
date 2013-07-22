find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})

add_library(flandmark_static STATIC flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
target_link_libraries(flandmark_static ${OpenCV_LIBS})
if(CMAKE_COMPILER_IS_GNUCC)
    set_target_properties(flandmark_static PROPERTIES COMPILE_FLAGS -fPIC)
endif(CMAKE_COMPILER_IS_GNUCC)
set_property(TARGET flandmark_static PROPERTY COMPILE_DEFINITIONS FLANDMARK_STATIC)

add_library(flandmark_shared SHARED flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
target_link_libraries(flandmark_shared ${OpenCV_LIBS})

#setup Config.cmake
SET(FLANDMARK_BASE_DIR "${PROJECT_SOURCE_DIR}/libflandmark")
set(FLANDMARK_BINARY_DIR "${PROJECT_BINARY_DIR}/libflandmark")
configure_file(flandmarkConfig.cmake.in
  "${PROJECT_BINARY_DIR}/libflandmark/flandmarkConfig.cmake" @ONLY)

export(PACKAGE flandmark)
