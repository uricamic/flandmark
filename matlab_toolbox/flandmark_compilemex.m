%% flandmark_compilemex.m
% Compilation script for flandmark mex-files.
% Tested with OpenCV version 2.3 (linux, )
% 
% 09-12-12, Michal Uricar

clc;
clearvars; close all;

% Change this to your version of opencv
OpenCVmajor = 2;
OpenCVminor = 4;

% Setup paths accpording to system platform
if (isunix&&~ismac)
    fprintf('Compiling mex-files on Unix platform...\n');
    
    include = ' -I/usr/include/opencv/ -I/usr/include/ ';
    libpath = '/usr/lib/';
    libflandmark = ' -L../bin/libflandmark/ -lflandmark_static ';
    
    if (OpenCVminor < 3)
        cxcore = dir([libpath 'libcxcore.so.' num2str(OpenCVmajor) '.' num2str(OpenCVminor)]);
        cxcore = cxcore.name;
        cv = dir([libpath 'libcv.so.' num2str(OpenCVmajor) '.' num2str(OpenCVminor)]);
        cv = cv.name;
        cvaux = dir([libpath 'libcvaux.so.' num2str(OpenCVmajor) '.' num2str(OpenCVminor)]);
        cvaux = cvaux.name;
        
        cvlibs = [libpath cxcore ' ' libpath cv ' ' libpath cvaux ' '];
    else
        opencv_core = dir([libpath 'libopencv_core.so.' num2str(OpenCVmajor) '.' num2str(OpenCVminor)]);
        opencv_core = opencv_core.name;
        opencv_imgproc = dir([libpath 'libopencv_imgproc.so.' num2str(OpenCVmajor) '.' num2str(OpenCVminor)]);
        opencv_imgproc = opencv_imgproc.name;
        
        cvlibs = [libpath opencv_core ' ' libpath opencv_imgproc ' '];
    end;
    
    % compile mexfiles
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_detector']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_base_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_detector_base']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_write_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_write_model']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_load_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_load_model']);
    eval(['mex -O -largeArrayDims ./mex/lbp2sparse_mex.cpp -I../libflandmark/' libflandmark ' -o ../matlab_toolbox/mex/lbp2sparse']);
    eval(['mex -O -largeArrayDims ./mex/argmax_mex.cpp -I../libflandmark/ ' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/argmax_mex']);
    eval(['mex -O -largeArrayDims ./mex/getNormalizedFrame_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/getNormalizedFrame']);

    fprintf('Compilation finished.\n');
end;

if (ispc)
    fprintf('Compiling mex-files on Windows platform...\n');
    
    include = [' -Ic:\OpenCV' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '\build\include\opencv\', ... 
               ' -Ic:\OpenCV' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '\build\include\' ];
    libpath = ['c:\OpenCV' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '\build\x64\vc10\lib\'];
    libflandmark = ' -L../bin7/libflandmark/Release -lflandmark_static ';
    
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
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_detector']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_base_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_detector_base']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_write_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_write_model']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_load_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/flandmark_load_model']);
    eval(['mex -O -largeArrayDims ./mex/lbp2sparse_mex.cpp -I../libflandmark/' libflandmark ' -output ../matlab_toolbox/mex/lbp2sparse']);
    eval(['mex -O -largeArrayDims ./mex/argmax_mex.cpp -I../libflandmark/ ' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/argmax_mex']);
    eval(['mex -O -largeArrayDims ./mex/getNormalizedFrame_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/getNormalizedFrame']);

    fprintf('Compilation finished.\n');
end;

if (ismac)
    fprintf('Compiling mex-files on Mac platform...\n');
    
    include = ' -I/opt/local/include/opencv/ -I/opt/local/include/ -I/usr/local/include/opencv/ -I/usr/local/include/'; 
    if(numel(dir('/usr/local/lib/libopencv*')))
        libpath = '/usr/local/lib/';
    else
        libpath = '/opt/local/lib/';
    end
    libflandmark = ' -L../bin/libflandmark/ -lflandmark_static ';
    
    % files = dir([libpath 'libopencv*.dylib']);  
    % e.g. libopencv_core.2.4.2.dylib
    
    if (OpenCVminor < 3)
        cxcore = dir([libpath 'libcxcore.' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '.dylib']);
        cxcore = cxcore.name;
        cv = dir([libpath 'libcv.' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '.dylib']);
        cv = cv.name;
        cvaux = dir([libpath 'libcvaux.' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '.dylib']);
        cvaux = cvaux.name;
        
        cvlibs = [libpath cxcore ' ' libpath cv ' ' libpath cvaux ' '];
    else
        opencv_core = dir([libpath 'libopencv_core.' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '.dylib']);
        opencv_core = opencv_core.name;
        opencv_imgproc = dir([libpath 'libopencv_imgproc.' num2str(OpenCVmajor) '.' num2str(OpenCVminor) '.dylib']);
        opencv_imgproc = opencv_imgproc.name;
        
        cvlibs = [libpath opencv_core ' ' libpath opencv_imgproc ' '];
    end;
    
    % compile mexfiles
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_detector']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_detector_base_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_detector_base']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_write_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_write_model']);
    eval(['mex -O -largeArrayDims ./mex/flandmark_load_model_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/flandmark_load_model']);
    eval(['mex -O -largeArrayDims ./mex/lbp2sparse_mex.cpp -I../libflandmark/' libflandmark ' -o ../matlab_toolbox/mex/lbp2sparse']);
    eval(['mex -O -largeArrayDims ./mex/argmax_mex.cpp -I../libflandmark/ ' include libflandmark cvlibs ' -output ../matlab_toolbox/mex/argmax_mex']);
    eval(['mex -O -largeArrayDims ./mex/getNormalizedFrame_mex.cpp -I../libflandmark/' include libflandmark cvlibs ' -o ../matlab_toolbox/mex/getNormalizedFrame']);

    fprintf('Compilation finished.\n');
end;
