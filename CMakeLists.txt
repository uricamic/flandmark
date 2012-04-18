CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(flandmark)

# The version number.
SET(flandmark_VERSION_MAJOR 1)
SET(flandmark_VERSION_MINOR 6)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

OPTION(BUILD_MATLAB_BINDINGS "Build Matlab bindings" ON)
SET(FLANDMARK_LEARNING_DIR "${PROJECT_SOURCE_DIR}/learning" CACHE STRING "Directory where will be copyied mex files")

ADD_SUBDIRECTORY(cpp)

# copy some important files to the binary directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/flandmark_model.dat ${CMAKE_CURRENT_BINARY_DIR}/cpp/flandmark_model.dat COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/face.jpg ${CMAKE_CURRENT_BINARY_DIR}/cpp/face.jpg COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/groupphoto.jpg ${CMAKE_CURRENT_BINARY_DIR}/cpp/groupphoto.jpg COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/seq_bruges04_300frames.avi ${CMAKE_CURRENT_BINARY_DIR}/cpp/seq_bruges04_300frames.avi COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/haarcascade_frontalface_alt.xml ${CMAKE_CURRENT_BINARY_DIR}/cpp/haarcascade_frontalface_alt.xml COPYONLY)