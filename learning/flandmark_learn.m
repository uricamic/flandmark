%% flandmark_learn.m
% TODO - add description
%
% 14-07-10 Michal Uricar
% 17-05-11 Michal Uricar
% 02-08-11 Michal Uricar, flandmark2 - detection of corners of eyes and mouth
% 21-03-12 Michal Uricar, LFW annotation in one file
% 13-09-12 Michal Uricar, refactoring

clc; 
close all; clearvars;

%% User defined options 

% setup paths to face database (and benchmarks)
ROOT_FACE_FOLDER = '../LFW/';                           % path to folder with database images
DB_ANN = '../data/LFW_annotation.txt';                  % annotation file for the database
%TRN_EXAMPLES = '../benchmarks/LFW_TRN_annotation.txt';  % benchmark db division, leave empty to create different division
%VAL_EXAMPLES = '../benchmarks/LFW_VAL_annotation.txt';  % benchmark db division, leave empty to create different division
%TST_EXAMPLES = '../benchmarks/LFW_TST_annotation.txt';  % benchmark db division, leave empty to create different division
TRN_EXAMPLES = [];
VAL_EXAMPLES = [];
TST_EXAMPLES = [];

% options configuration
% order of components in LFW database
options.compnames = { 'center', 'canthus-rl', 'canthus-lr', 'mouth-corner-r', 'mouth-corner-l', ...
                      'canthus-rr', 'canthus-ll', 'nose' };
options.bw = [40; 40];
options.components = [20 10 10 10 10 10 10 10;
                      20 10 10 10 10 10 10 10];
options.comselect = [10, 4, 5, 8, 9, 3, 6, 10];
% number of components
options.M = size(options.components, 2);
% base window margin (percents)
options.bw_margin = [20; 20];
% path to LFW database (where images and anotations are stored) must contain subdir /mgt_v2 with xml annotations
options.image_path = ROOT_FACE_FOLDER;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add paths

	addpath('./code/');

%% Check project dir structure

    fprintf('Checking project dir structure...\n');
    
    DIRS = {'./code/img/', './code/results/', './code/MAT'};
    
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
    
    DIRS = {'./bmrm', './libqp', './code/Functions', ...
            ROOT_FACE_FOLDER};
    
    for i = 1 : numel(DIRS)
        fprintf('%s ... ', DIRS{i});
        if (~exist(DIRS{i}, 'dir'))
            fprintf('MISSING\n');
            error('??? Corrupted project structure.');
        else 
            fprintf('OK\n');
        end;
    end;
    
%% Save above specified user data

    options.image_path = ['../' ROOT_FACE_FOLDER];

    save('./code/MAT/paths.mat', 'ROOT_FACE_FOLDER', 'TRN_EXAMPLES', 'VAL_EXAMPLES', 'TST_EXAMPLES', 'DB_ANN');
    save('./code/MAT/options.mat', 'options');
    
%% Compile mex files if neccessary

    fprintf('Checking for neccessary mex-files...\n');

    if (~exist(['../matlab_toolbox/mex/lbp2sparse.' mexext], 'file') || ...
        ~exist(['../matlab_toolbox/mex/argmax_mex.' mexext], 'file') || ...
        ~exist(['../matlab_toolbox/mex/getNormalizedFrame.' mexext], 'file') || ...
        ~exist(['../matlab_toolbox/mex/flandmark_detector.' mexext], 'file'))
        error(['??? Mex-files not found! Please compile the flandmark with BUILD_MATLAB_BINDINGS turned on'...
               ' or run /matlab_toolbox/flandmark_compilemex.m to create them.\n']);
    end
    
    fprintf('Done.\n');
    
%% Collect statistics about components positions

    % check if statistics already exist
    if (exist('./code/results/statistics_orig.mat', 'file'))
        fprintf('Statistics found. Skipping...\n');
    else 
        fprintf('Collecting statistics about components...\n');
        
        run './code/runstats.m';
    end;

%% Division of database into 3 parts (60% traning, 20% validating, 20% testing)

    % check if division already done.
    if (exist('./code/MAT/TRN.mat', 'file') && ...
            exist('./code/MAT/VAL.mat', 'file') && ...
            exist('./code/MAT/TST.mat', 'file'))
        fprintf('Division already done. Skipping...\n');
    else
        fprintf('Dividing database into 3 parts...\n');
        
        run './code/dbDivision.m';
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
    save('./code/MAT/functions_gdisp.mat', 'f_argmax', 'f_risk', ...
        'f_buildPsi', 'f_createMapTable', 'displacement');

    %% Configuration - nImages, height_of_pyramid
    fprintf('Looking for config file...\n');
    if (~exist('./code/MAT/config.mat', 'file'))
        fprintf('Not found. Creating...\nconfig = \n');
        config.nImages = 'all';
        config.heightOfPyramid = 4;
        disp(config);
        save('./code/MAT/config.mat', 'config');
    else 
        fprintf('File %s found.\n', './code/MAT/config.mat');
    end;

    fprintf('Prepare data for LBP\n');
    run './code/prepDataForLBP.m';
    
%% Call BMRM solver

    % prepare lamdas
    fprintf('Looking for file with lambda values...\n');
    if (~exist('./code/MAT/lambdas.mat', 'file'))
        lambdas = [1e-3, 1e-2, 1e-1, 1];
        fprintf('Not found. Creating...\n');
        save('./code/MAT/lambdas.mat', 'lambdas');
        disp(lambdas);
    else
        fprintf('Lambda values found:');
        load('./code/MAT/lambdas.mat');
        disp(lambdas);
    end;

    % prepare normalization coefficients for each image in database
    if (~exist('./code/MAT/kappas.mat', 'file'))
	run('./code/computeKappas.m');
    end;
    
    fprintf('Calling BMRM solver...\n');
    run './code/test_bmrm.m';
    
%% Validation

    if (exist('./code/results/W.mat', 'file'))
        fprintf('Best vector W already found. Skipping...\n');
    else
        fprintf('Running validaton...\n');
        
        run './code/testValidation.m';
    end;
    
%% Test

    if (exist('./code/results/errors.mat', 'file'))
        fprintf('Test already passed. Skipping...\n');
    else
        fprintf('Running test...\n');
        
        run './code/test1.m';
    end;
    
    fprintf('Show histograms.\n');
    run './code/test2.m';
    
%% Generate model.mat (for use with flandmark_load_model - binary dat file creation)

    run './code/detectorModelBuild.m';
    
%% Timestamp
      
 fprintf(1,'Finished on %s.\n', datestr(now));
 