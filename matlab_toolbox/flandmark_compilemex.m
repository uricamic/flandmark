%% flandmark_compilemex.m
% Compilation script for flandmark mex-files.
% Tested with OpenCV version 2.3 (linux, )
% 
% 09-12-12, Michal Uricar

clc;
clearvars; close all;

% Change this to your version of opencv
OpenCVmajor = 24;
OpenCVminor = 10;

% Setup paths accpording to system platform
    fprintf('Compiling mex-files on Windows platform...\n');
    include = [' -Ie:\opencv\build\include\opencv\' ' -Ie:\opencv\build\include\' ];
    libpath = 'e:\opencv\build\x64\vc10\lib\';
    libflandmark = ' -L../flandmark/libflandmark_static/Release -lflandmark_static ';
    
    if (OpenCVminor < 3)
        cxcore = dir([libpath 'cxcore' num2str(OpenCVmajor) num2str(OpenCVminor) '*.lib']);
        cxcore = cxcore.name;
        cv = dir([libpath 'cv' num2str(OpenCVmajor) num2str(OpenCVminor) '*.lib']);
        cv = cv.name;
        cvaux = dir([libpath 'cvaux' num2str(OpenCVmajor) num2str(OpenCVminor) '*.lib']);
        cvaux = cvaux.name;
        
        cvlibs = [libpath cxcore ' ' libpath cv ' ' libpath cvaux ' '];
    else
        opencv_core = dir([libpath 'opencv_core' num2str(OpenCVmajor) num2str(OpenCVminor) '*.lib']);
        opencv_core = opencv_core.name;
        opencv_imgproc = dir([libpath 'opencv_imgproc' num2str(OpenCVmajor) num2str(OpenCVminor) '*.lib']);
        opencv_imgproc = opencv_imgproc.name;
        
        cvlibs = [libpath opencv_core ' ' libpath opencv_imgproc ' '];
    end;
    
    % compile mexfiles
%     ['mex -O -largeArrayDims ./mex/flandmark_detector_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_detector']
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_detector']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_base_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_detector_base']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_write_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_write_model']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_load_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_load_model']);
    eval(['mex -O -largeArrayDims ./mex/lbp2sparse_mex.cpp -I../libflandmark/' libflandmark ' -output ../matlab_toolbox/mex/lbp2sparse']);
    eval(['mex -O -largeArrayDims ./mex/argmax_mex.cpp -I../libflandmark/ ' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/argmax_mex']);
    eval(['mex -O -largeArrayDims ./mex/getNormalizedFrame_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/getNormalizedFrame']);

    fprintf('Compilation finished.\n');

