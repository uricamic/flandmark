find_package( OpenCV REQUIRED )
include_directories(${OPENCV_INCLUDE_DIR})

if(WIN32)
	find_program(MEX_CMD mex.bat)
	find_program(MEXEXT_CMD mexext.bat)
else(WIN32)
	find_program(MEX_CMD mex)
	find_program(MEXEXT_CMD mexext)
endif(WIN32)

set(MEX_NAME_flandmark_detector flandmark_detector)
set(MEX_NAME_flandmark_detector_base flandmark_detector_base)
set(MEX_NAME_flandmark_write_model flandmark_write_model)
set(MEX_NAME_flandmark_load_model flandmark_load_model)
set(MEX_NAME_lbp2sparse lbp2sparse)
set(MEX_NAME_argmax argmax_mex)
set(MEX_NAME_getNF getNormalizedFrame)

if(${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 -lcxcore -lcv -lcvaux)
	elseif (WIN32)
		set(CV_LIBS_1 -lcxcore${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} -lcv${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} -lcvaux${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
else(${OpenCV_VERSION_MINOR} LESS 3)
	if (UNIX)
		set(CV_LIBS_1 -lopencv_core -lopencv_imgproc)
	elseif (WIN32)
		set(CV_LIBS_1 -lopencv_core${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH} -lopencv_imgproc${OpenCV_VERSION_MAJOR}${OpenCV_VERSION_MINOR}${OpenCV_VERSION_PATCH})
	elseif (APPLE)
	endif (UNIX)
endif(${OpenCV_VERSION_MINOR} LESS 3)

if(MEX_CMD AND MEXEXT_CMD)
	get_filename_component(MEX_REAL_CMD ${MEX_CMD} ABSOLUTE)
	get_filename_component(MEX_PATH ${MEX_REAL_CMD} PATH)

	get_filename_component(MEXEXT_REAL_CMD ${MEXEXT_CMD} ABSOLUTE)
	get_filename_component(MEXEXT_PATH ${MEXEXT_REAL_CMD} PATH)

	#message(STATUS "Found MATLAB at: " ${MEX_PATH})
	#message(STATUS "OpenCV_INSTALL_PATH =" ${OpenCV_INSTALL_PATH})

	if(MEX_PATH STREQUAL MEXEXT_PATH)

		execute_process(COMMAND ${MEXEXT_REAL_CMD} OUTPUT_VARIABLE MEX_EXTENSION OUTPUT_STRIP_TRAILING_WHITESPACE)
		set(MEX_FILE_flandmark_detector ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_detector}.${MEX_EXTENSION})
		set(MEX_FILE_flandmark_detector_base ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_detector_base}.${MEX_EXTENSION})
		set(MEX_FILE_flandmark_write_model ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_write_model}.${MEX_EXTENSION})
		set(MEX_FILE_flandmark_load_model ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_load_model}.${MEX_EXTENSION})
		set(MEX_FILE_lbp2sparse ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_lbp2sparse}.${MEX_EXTENSION})
		set(MEX_FILE_argmax ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_argmax}.${MEX_EXTENSION})
		set(MEX_FILE_getNF ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_getNF}.${MEX_EXTENSION})

		find_path(OpenCV_INCLUDE "cv.h" PATHS "${OpenCV_DIR}" PATH_SUFFIXES "include" "include/opencv" DOC "Include directory")
		#message("-- cflags opencv: -I" ${OpenCV_INCLUDE})
		#message("-- OpenCV_DIR = ${OpenCV_DIR}")

		if (WIN32)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_detector}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_flandmark_detector}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_detector_base}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_base_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_flandmark_detector_base}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector_base}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector_base} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_write_model}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_write_model_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_flandmark_write_model}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_write_model}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_write_model} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_load_model}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_load_model_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_flandmark_load_model}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_load_model}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_load_model} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_lbp2sparse}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/lbp2sparse_mex.cpp ${PROJECT_SOURCE_DIR}/libflandmark/liblbp.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -output ${MEX_NAME_lbp2sparse}
				COMMENT "Building MEX extension ${MEX_FILE_lbp2sparse} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_argmax}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/argmax_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_argmax}
				COMMENT "Building MEX extension ${MEX_FILE_argmax} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_getNF}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/getNormalizedFrame_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -output ${MEX_NAME_getNF}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_getNF}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_getNF} "
			)

		else (WIN32)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_detector}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_flandmark_detector}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_detector_base}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_base_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_flandmark_detector_base}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector_base}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector_base} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_write_model}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_write_model_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_flandmark_write_model}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_write_model}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_write_model} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_flandmark_load_model}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_load_model_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_flandmark_load_model}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_load_model}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_flandmark_load_model} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_lbp2sparse}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/lbp2sparse_mex.cpp ${PROJECT_SOURCE_DIR}/libflandmark/liblbp.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -o ${MEX_NAME_lbp2sparse}
				COMMENT "Building MEX extension ${MEX_FILE_lbp2sparse} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_argmax}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/argmax_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_argmax}
				COMMENT "Building MEX extension ${MEX_FILE_argmax} "
			)
			ADD_CUSTOM_COMMAND(
				OUTPUT ${MEX_FILE_getNF}
				COMMAND ${MEX_REAL_CMD}
				ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/getNormalizedFrame_mex.cpp -I${PROJECT_SOURCE_DIR}/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static ${CV_LIBS_1} -o ${MEX_NAME_getNF}
				DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_getNF}_mex.cpp
				COMMENT "Building MEX extension ${MEX_FILE_getNF} "
			)
		endif (WIN32)

		add_custom_target(${MEX_NAME_flandmark_detector} ALL DEPENDS ${MEX_FILE_flandmark_detector})
		add_custom_target(${MEX_NAME_flandmark_detector_base} ALL DEPENDS ${MEX_FILE_flandmark_detector_base})
		add_custom_target(${MEX_NAME_flandmark_write_model} ALL DEPENDS ${MEX_FILE_flandmark_write_model})
		add_custom_target(${MEX_NAME_flandmark_load_model} ALL DEPENDS ${MEX_FILE_flandmark_load_model})
		add_custom_target(${MEX_NAME_lbp2sparse} ALL DEPENDS ${MEX_FILE_lbp2sparse})
		add_custom_target(${MEX_NAME_argmax} ALL DEPENDS ${MEX_FILE_argmax})
		add_custom_target(${MEX_NAME_getNF} ALL DEPENDS ${MEX_FILE_getNF})

		# copy mex-files
		add_custom_command(
			TARGET ${MEX_NAME_getNF}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_getNF} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_getNF}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_argmax}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_argmax} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_argmax}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_lbp2sparse}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_lbp2sparse} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_lbp2sparse}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_flandmark_detector}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_detector} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_flandmark_detector}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_flandmark_detector_base}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_detector_base} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_flandmark_detector_base}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_flandmark_write_model}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_write_model} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_flandmark_write_model}"
			)
		add_custom_command(
			TARGET ${MEX_NAME_flandmark_load_model}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_load_model} ${PROJECT_SOURCE_DIR}/matlab_toolbox/mex
			COMMENT "Copying ${MEX_FILE_flandmark_load_model}"
			)

	else()
		message(WARNING "The 'mex' and 'mexext' programs have been found in different locations. It's likely that one of them is not part of the MATLAB installation. Make sure that the 'bin' directory from the MATLAB installation is in PATH")
		set(BUILD_MATLAB_BINDINGS OFF)
	endif()
endif()
