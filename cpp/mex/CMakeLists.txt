FIND_PACKAGE( OpenCV REQUIRED )
INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIR})

if(WIN32)
    find_program(MEX_CMD mex.bat)
    find_program(MEXEXT_CMD mexext.bat)
else()
    find_program(MEX_CMD mex)
    find_program(MEXEXT_CMD mexext)
endif()

SET(MEX_NAME_flandmark_detector flandmark_detector)
SET(MEX_NAME_flandmark_detector_base flandmark_detector_base)
SET(MEX_NAME_flandmark_load_model flandmark_load_model)
SET(MEX_NAME_lbp2sparse lbp2sparse)
SET(MEX_NAME_argmax argmax_mex)
SET(MEX_NAME_getNF getNormalizedFrame)

if(MEX_CMD AND MEXEXT_CMD)
    get_filename_component(MEX_REAL_CMD ${MEX_CMD} ABSOLUTE)
    get_filename_component(MEX_PATH ${MEX_REAL_CMD} PATH)

    get_filename_component(MEXEXT_REAL_CMD ${MEXEXT_CMD} ABSOLUTE)
    get_filename_component(MEXEXT_PATH ${MEXEXT_REAL_CMD} PATH)

    message(STATUS "Found MATLAB at: " ${MEX_PATH})
    message(STATUS "OpenCV_INSTALL_PATH =" ${OpenCV_INSTALL_PATH})

    if(MEX_PATH STREQUAL MEXEXT_PATH)

        EXECUTE_PROCESS(COMMAND ${MEXEXT_REAL_CMD} OUTPUT_VARIABLE MEX_EXTENSION OUTPUT_STRIP_TRAILING_WHITESPACE)
        SET(MEX_FILE_flandmark_detector ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_detector}.${MEX_EXTENSION})
        SET(MEX_FILE_flandmark_detector_base ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_detector_base}.${MEX_EXTENSION})
        SET(MEX_FILE_flandmark_load_model ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_flandmark_load_model}.${MEX_EXTENSION})
        SET(MEX_FILE_lbp2sparse ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_lbp2sparse}.${MEX_EXTENSION})
        SET(MEX_FILE_argmax ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_argmax}.${MEX_EXTENSION})
        SET(MEX_FILE_getNF ${CMAKE_CURRENT_BINARY_DIR}/${MEX_NAME_getNF}.${MEX_EXTENSION})

	if(${OpenCV_VERSION_MINOR} LESS 3)
       
                find_path(OpenCV_INCLUDE "cv.h" PATHS "${OpenCV_DIR}"
                        PATH_SUFFIXES "include/opencv" DOC "Include
                        directory")

                message("--cflags opencv: -I" ${OpenCV_INCLUDE})
                
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_detector}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims
	            ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_mex.cpp
	            -I${PROJECT_SOURCE_DIR}/cpp/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lcxcore -lcv -lcvaux -o ${MEX_NAME_flandmark_detector}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_detector_base}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims
	            ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_base_mex.cpp
	            -I${PROJECT_SOURCE_DIR}/cpp/libflandmark -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lcxcore -lcv -lcvaux -o ${MEX_NAME_flandmark_detector_base}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector_base}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector_base}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_load_model}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims
	            ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_load_model_mex.cpp
	            -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lcxcore -lcv -lcvaux -o ${MEX_NAME_flandmark_load_model}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_load_model}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_load_model}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_lbp2sparse}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/lbp2sparse_mex.cpp ${PROJECT_SOURCE_DIR}/cpp/libflandmark/liblbp.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -o ${MEX_NAME_lbp2sparse}
	            COMMENT "Building MEX extension ${MEX_FILE_lbp2sparse}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_argmax}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/argmax_mex.cpp -o ${MEX_NAME_argmax}
	            COMMENT "Building MEX extension ${MEX_FILE_argmax}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_getNF}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims
	            ${CMAKE_CURRENT_SOURCE_DIR}/getNormalizedFrame_mex.cpp
	            -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -I${OpenCV_INCLUDE} -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lcxcore -lcv -lcvaux -o ${MEX_NAME_getNF}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_getNF}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_getNF}"
	        )

        else(${OpenCV_VERSION_MINOR} LESS 3)

	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_detector}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_mex.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark -I${OpenCV_INSTALL_PATH}/include/opencv -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lopencv_core -lopencv_imgproc -o ${MEX_NAME_flandmark_detector}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_detector_base}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_detector_base_mex.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark -I${OpenCV_INSTALL_PATH}/include/opencv -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lopencv_core -lopencv_imgproc -o ${MEX_NAME_flandmark_detector_base}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_detector_base}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_detector_base}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_flandmark_load_model}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/flandmark_load_model_mex.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -I${OpenCV_INSTALL_PATH}/include/opencv -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lopencv_core -lopencv_imgproc -o ${MEX_NAME_flandmark_load_model}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_flandmark_load_model}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_flandmark_load_model}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_lbp2sparse}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/lbp2sparse_mex.cpp ${PROJECT_SOURCE_DIR}/cpp/libflandmark/liblbp.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -o ${MEX_NAME_lbp2sparse}
	            COMMENT "Building MEX extension ${MEX_FILE_lbp2sparse}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_argmax}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/argmax_mex.cpp -o ${MEX_NAME_argmax}
	            COMMENT "Building MEX extension ${MEX_FILE_argmax}"
	        )
	
	        ADD_CUSTOM_COMMAND(
	            OUTPUT ${MEX_FILE_getNF}
	            COMMAND ${MEX_REAL_CMD}
	            ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/getNormalizedFrame_mex.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -I${OpenCV_INSTALL_PATH}/include/opencv -L${PROJECT_BINARY_DIR}/cpp/libflandmark -lflandmark_static -lopencv_core -lopencv_imgproc -o ${MEX_NAME_getNF}
	            #ARGS -largeArrayDims ${CMAKE_CURRENT_SOURCE_DIR}/getNormalizedFrame_mex.cpp -I${PROJECT_SOURCE_DIR}/cpp/libflandmark/ -I${OpenCV_INSTALL_PATH}/include/opencv -L${PROJECT_BINARY_DIR}/libflandmark -lflandmark_static -lopencv_core -lopencv_imgproc -lopencv_highgui -o ${MEX_NAME_getNF}
	            DEPENDS flandmark_static ${CMAKE_CURRENT_SOURCE_DIR}/${MEX_NAME_getNF}_mex.cpp
	            COMMENT "Building MEX extension ${MEX_FILE_getNF}"
	        )

	endif(${OpenCV_VERSION_MINOR} LESS 3)

        add_custom_target(${MEX_NAME_flandmark_detector} ALL DEPENDS ${MEX_FILE_flandmark_detector})
        add_custom_target(${MEX_NAME_flandmark_detector_base} ALL DEPENDS ${MEX_FILE_flandmark_detector_base})
        add_custom_target(${MEX_NAME_flandmark_load_model} ALL DEPENDS ${MEX_FILE_flandmark_load_model})
        add_custom_target(${MEX_NAME_lbp2sparse} ALL DEPENDS ${MEX_FILE_lbp2sparse})
        add_custom_target(${MEX_NAME_argmax} ALL DEPENDS ${MEX_FILE_argmax})
        add_custom_target(${MEX_NAME_getNF} ALL DEPENDS ${MEX_FILE_getNF})

#		macro(copy_files GLOBPAT DESTINATION)
##file(GLOB COPY_FILES RELATIVE ${CMAKE_CURRENT_BINARY_DIR} ${GLOBPAT})
#			file(GLOB COPY_FILES RELATIVE ${CMAKE_BINARY_DIR}/cpp/mex ${GLOBPAT})
#			add_custom_target(
#	    		copy ALL
#	    		COMMENT "Copying files: ${GLOBPAT}: ${COPY_FILES}"
#			)
#			foreach(FILENAME ${COPY_FILES})
##set(SRC "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
#				set(SRC "${CMAKE_BINARY_DIR}/cpp/mex/${FILENAME}")
#	  			set(DST "${DESTINATION}/${FILENAME}")
#	  			add_custom_command(
#	    			TARGET copy
#	    			COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
#					COMMENT "Copying file ${SRC} to ${DST}."
#	    		)
#			endforeach(FILENAME)
#   	endmacro(copy_files)
#		copy_files("${CMAKE_BINARY_DIR}/cpp/mex/*.${MEX_EXTENSION}" "${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex")


#	macro(copy_files GLOBPAT DESTINATION)
#		file(GLOB COPY_FILES RELATIVE ${CMAKE_CURRENT_BINARY_DIR} ${GLOBPAT})
#		add_custom_target(
#			copy ALL
#			COMMENT "Copying files: ${GLOBPAT}"
#			COMMENT "COPY_FILES: ${COPY_FILES}")
#		foreach(FILENAME ${COPY_FILES})
#			set(SRC "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
#			set(DST "${DESTINATION}/${FILENAME}")
#			add_custom_command(
#				TARGET copy
#				COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
#				)
#		endforeach(FILENAME)
#	endmacro(copy_files)

#	copy_files("${CMAKE_CURRENT_BINARY_DIR}/*.${MEX_EXTENSION}" "${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex")


	# copy mex-files
	add_custom_command(
			TARGET ${MEX_NAME_getNF}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_getNF} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_getNF}"
	)
	add_custom_command(
			TARGET ${MEX_NAME_argmax}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_argmax} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_argmax}"
	)
	add_custom_command(
			TARGET ${MEX_NAME_lbp2sparse}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_lbp2sparse} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_lbp2sparse}"
	)
	add_custom_command(
			TARGET ${MEX_NAME_flandmark_detector}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_detector} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_flandmark_detector}"
	)
	add_custom_command(
			TARGET ${MEX_NAME_flandmark_detector_base}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_detector_base} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_flandmark_detector_base}"
	)
	add_custom_command(
			TARGET ${MEX_NAME_flandmark_load_model}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy ${MEX_FILE_flandmark_load_model} ${PROJECT_SOURCE_DIR}/learning/flandmark_data/libocas/mex
			COMMENT "Copying ${MEX_FILE_flandmark_load_model}"
	)

    else()
        message(WARNING "The 'mex' and 'mexext' programs have been found in different locations. It's likely that one of them is not part of the MATLAB instalation. Make sure that the 'bin' directory from the MATLAB instalation is in PATH")
        set(BUILD_MATLAB_BINDINGS OFF)
    endif()
endif()
