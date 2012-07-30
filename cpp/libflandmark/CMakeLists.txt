FIND_PACKAGE( OpenCV REQUIRED )
INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIR})

if(${OpenCV_VERSION_MINOR} LESS 3)
        
    add_library(flandmark_static STATIC flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
    target_link_libraries(flandmark_static cxcore cv cvaux)
    if(CMAKE_COMPILER_IS_GNUCC)
        set_target_properties(flandmark_static PROPERTIES COMPILE_FLAGS -fPIC)
    endif()
    set_property(TARGET flandmark_static PROPERTY COMPILE_DEFINITIONS FLANDMARK_STATIC)

    ADD_LIBRARY(flandmark_shared SHARED flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
    target_link_libraries(flandmark_shared cxcore cv cvaux)

else(${OpenCV_VERSION_MINOR} LESS 3)

    add_library(flandmark_static STATIC flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
    target_link_libraries(flandmark_static opencv_core opencv_imgproc)
    if(CMAKE_COMPILER_IS_GNUCC)
        set_target_properties(flandmark_static PROPERTIES COMPILE_FLAGS -fPIC)
    endif()
    set_property(TARGET flandmark_static PROPERTY COMPILE_DEFINITIONS FLANDMARK_STATIC)
    
    ADD_LIBRARY(flandmark_shared SHARED flandmark_detector.cpp flandmark_detector.h liblbp.cpp liblbp.h)
    target_link_libraries(flandmark_shared opencv_core opencv_imgproc)

endif(${OpenCV_VERSION_MINOR} LESS 3)
