% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% flandmark_learn.m
% 
% Main script for learning parameters of facial landmark detector based on
% Deformable Part Models (DPM) learned by structured output SVM algorithm
% 
% Script assumes following dir structure:
% ./flandmark_data/Functions/   - directory containing helper functinos
%                                 (such as getImageFrame, etc.)
% ./flandmark_data/MAT/         - directory with neccessary .mat files
% ./flandmark_data/results/     - directory for storing results of some
%                                 tests 
% ./flandmark_data/LFW/         - symbolic link to LFW database directory
%                                 (must also contain /mgt/ subdir) 
% 
% Directory structure is checked by the script.
% 
% First step is pruning of original image database to get rid of incomplete
% or corrupted data. Incomplete are images with incomplete XML annotation
% (for example one landmark is missing) or images without detected bbox of
% face, corrupted are images with multiple faces detected or with multiple
% faces in one XML annotation.
% 
% Second step is collecting statistics of components positions. This means
% that every image in pruned database is cropped to detected bbox extended
% by parameter options.bw_margin and than rescaled to size of parameter
% option.bw (base window). This we call normalized frame. Each component
% has assigned its own bbox (options.components). We count images that has
% all components in normalized frame. Parameters are set to give good
% number of inliers (for current parameters about 97.7% of images passed).
% Result of this test is matrix S holding areas of possible positions of
% each landmark point (component)
% 
% Next step is separation of face database into 3 parts (training,
% validating and testing)
% 
% With properly selected zero component next step prepares data for
% classifier learning--- in this step data for computing the 
% LBP features on all 3 parts of database (TRN, VAL, TST) are computed.
% 
% Next step is classifier learning done by BMRM solver for some lambda
% values.
% 
% Best vector W is selected in process of validation--- it is the one with
% the smallest mean error.
% 
% Last step is collection of test statistics. 
%
% 14-07-10 Michal Uricar
% 17-05-11 Michal Uricar
% 02-08-11 Michal Uricar, flandmark2 - detection of corners of eyes and mouth
% 21-03-12 Michal Uricar, LFW annotation in one file

clc; 
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add paths

	addpath('./flandmark_data/');

%% Check project dir structure

    fprintf('Checking project dir structure...\n');
    
    DIRS = {'./flandmark_data/img/', './flandmark_data/results/'};
    
    for i = 1 : numel(DIRS)
        fprintf('%s ... ', DIRS{i});
        if (~exist(DIRS{i}, 'dir'))
            fprintf('MISSING\n');
            fprintf('Creating dir %s...\n', DIRS{i});
            mkdir(DIRS{i});
        else
            fprintf('OK\n');
        end;
    end;
    
    DIRS = {'./flandmark_data/libocas', './flandmark_data/bmrm', './flandmark_data/Functions', ...
        './flandmark_data/MAT', './flandmark_data/LFW'};
    
    for i = 1 : numel(DIRS)
        fprintf('%s ... ', DIRS{i});
        if (~exist(DIRS{i}, 'dir'))
            fprintf('MISSING\n');
            error('??? Corrupted project structure.');
        else 
            fprintf('OK\n');
        end;
    end;
    
%% Compile mex files if neccessary

    fprintf('Checking for neccessary mex-files...\n');

    if (~exist(['./flandmark_data/libocas/mex/lbp2sparse.' mexext], 'file') || ...
        ~exist(['./flandmark_data/libocas/mex/argmax_mex.' mexext], 'file') || ...
        ~exist(['./flandmark_data/libocas/mex/getNormalizedFrame.' mexext], 'file') || ...
        ~exist(['./flandmark_data/libocas/mex/flandmark_detector.' mexext], 'file'))
        error('??? Mex-files not found! Please compile the flandmark with BUILD_MATLAB_BINDINGS turned on.\n');
    end
    
    fprintf('Done.\n');
    
%% Collect statistics about components positions

    % check if statistics already exist
    if (exist('./flandmark_data/results/statistics_orig.mat', 'file'))
        fprintf('Statistics found. Skipping...\n');
    else 
        fprintf('Collecting statistics about components...\n');
        
        run './flandmark_data/runstats.m';
    end;

%% Division of database into 3 parts (60% traning, 20% validating, 20% testing)

    % check if division already done.
    if (exist('./flandmark_data/MAT/TRN.mat', 'file') && ...
            exist('./flandmark_data/MAT/VAL.mat', 'file') && ...
            exist('./flandmark_data/MAT/TST.mat', 'file'))
        fprintf('Division already done. Skipping...\n');
    else
        fprintf('Dividing database into 3 parts...\n');
        
        run './flandmark_data/dbDivision.m';
    end;
    
%% Create data for BMRM solver (TRN, VAL, TST)

    %% Functions selection
    
    fprintf('Prepare functions configurations\n');
    f_argmax = @argmax2;
    f_risk = @risk2;
%     f_risk = @risk2_unnormalized;
    f_buildPsi = @buildPsi2;
    f_createMapTable = @createMapTable2;
    displacement = true;
    save('./flandmark_data/MAT/functions_gdisp.mat', 'f_argmax', 'f_risk', ...
        'f_buildPsi', 'f_createMapTable', 'displacement');

    %% Configuration - nImages, height_of_pyramid
    fprintf('Looking for config file...\n');
    if (~exist('./flandmark_data/MAT/config.mat', 'file'))
        fprintf('Not found. Creating...\nconfig = \n');
        config.nImages = 'all';
        config.heightOfPyramid = 4;
        disp(config);
        save('./flandmark_data/MAT/config.mat', 'config');
    else 
        fprintf('File %s found.\n', './flandmark_data/MAT/config.mat');
    end;

    fprintf('Prepare data for LBP\n');
    run './flandmark_data/prepDataForLBP.m';
    
%% Call BMRM solver

    % prepare lamdas
    fprintf('Looking for file with lambda values...\n');
    if (~exist('./flandmark_data/MAT/lambdas.mat', 'file'))
        lambdas = [1e-3, 1e-2, 1e-1, 1];
        fprintf('Not found. Creating...\n');
        save('./flandmark_data/MAT/lambdas.mat', 'lambdas');
        disp(lambdas);
    else
        fprintf('Lambda values found:');
        load('./flandmark_data/MAT/lambdas.mat');
        disp(lambdas);
    end;

    % prepare normalization coefficients for each image in database
    if (~exist('./flandmark_data/MAT/kappas.mat', 'file'))
	run('./flandmark_data/computeKappas.m');
    end;
    
    fprintf('Calling BMRM solver...\n');
    run './flandmark_data/test_bmrm.m';
    
%% Validation

    if (exist('./flandmark_data/results/W.mat', 'file'))
        fprintf('Best vector W already found. Skipping...\n');
    else
        fprintf('Running validaton...\n');
        
        run './flandmark_data/testValidation.m';
    end;
    
%% Test

    if (exist('./flandmark_data/results/errors.mat', 'file'))
        fprintf('Test already passed. Skipping...\n');
    else
        fprintf('Running test...\n');
        
        run './flandmark_data/test1.m';
    end;
    
    fprintf('Show histograms.\n');
    run './flandmark_data/test2.m';
    
%% Generate model.mat (for use with flandmark_load_model - binary dat file creation)

    run './flandmark_data/detectorModelBuild.m';
    
%% Timestamp
      
 fprintf(1,'Finished on %s.\n', datestr(now));
 