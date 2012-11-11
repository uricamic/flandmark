%% test_bmrm.m
% This script calls BMRM solver for given values of lambda on training
% (TRN) data. Output are corresponding vectors W_{\lambda}.
% 
% 30-07-10 Michal Uricar 

clc; 
clearvars; close all;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');
addpath('../libqp/matlab/');
addpath('../bmrm/');
addpath('../../matlab_toolbox/mex/');

%% Select functions

load('./MAT/functions_gdisp.mat');

%% Load data

fprintf('Loading data...\n');

% load config
load('./MAT/config.mat');

lbpconfig_fname = ['./MAT/lbpconfig_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
images_fname = ['./MAT/Images_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GT_nImages-' config.nImages '.mat'];

tic
load(lbpconfig_fname);
load(images_fname);
load(gt_fname);
load('./results/options.mat');
S = options.S;

switch (config.nImages)
    case 'all'
        nImages = size(Images, 2);
    case 'half'
        nImages = size(Images, 2)/2;
    otherwise
        nImages = str2num(config.nImages);
end;

load('./MAT/kappas.mat');

fprintf('Data loaded succesfully.\n');
toc

%% build structure data

fprintf('Building structure data...\n');

tic
data.lbp = lbpconfig; 
data.imSize = imSize; 
data.Images = Images; 
data.nImages = nImages;
options.S = S; data.Y = gt; data.options = options;
data.tmpData = getPsiMat(data, 1);
if (displacement)
    [data.options.PsiGS0, data.options.PsiGS1, data.options.PsiGS2] = getDisplacements(data.options);
end;
clear lbpconfig imSize Images nImages gt gt_fname images_fname lbpconfig_fname options S

% create mapTable for extracting components of W
data.mapTable = f_createMapTable(data.options, data.tmpData);
data.kappa = kappa;

fprintf('Structure data built.\n');
toc

%% BMRM

% lambda
% lambdas = [1e-2, 1e-1, 1];
load('./MAT/lambdas.mat');
lambdas = fliplr(lambdas)

% verbose on
options.verb = 1;
options.TolRel = 1e-2;

% fprintf('Starting MATLAB Pool...\n');
% try
%    matlabpool;
% catch ME1
%    idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Already connected to matlab pool. Skipping...\n');
%    end;
% end

for lambda = lambdas
    % check for W
    fname = ['./results/W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '_lambda_' num2str(lambda) '.mat'];
    if (exist(fname, 'file'))
        fprintf('Vector W found (file %s). Skipping...\n', fname);
    else
        fprintf('Lambda = %f. Call BMRM solver and save result to file %s.\n', lambda, fname);
        % call the BMRM solver
        [W, stat]= bmrm(data, f_risk, lambda, options);
        % save result
        save(fname, 'W', 'stat');
    end;
end;

% try
%    matlabpool close;
% catch ME2
%    idSegLast = regexp(ME2.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Cannot close matlabpool. Skipping...\n');
%    end;
% end;

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));
