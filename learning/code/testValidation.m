%% testValidation.m
% This script does validation of vectors W acquired from training data
% (see test_bmrm.m)
% 
% 30-07-10 Michal Uricar

% clc; 
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Select functions

load('./MAT/functions_gdisp.mat');

%% Load VAL data

% load config
load('./MAT/config.mat');
lbpconfig_fname = ['./MAT/lbpconfig_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
images_fname = ['./MAT/ImagesVal_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GTVal_nImages-' config.nImages '.mat'];

fprintf('Loading VAL data...\n');
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

tic
data.lbp = lbpconfig; data.imSize = imSize; data.Images = Images; data.nImages = nImages;
options.S = S; data.Y = gt; data.options = options;
data.tmpData = getPsiMat(data, 1);
if (displacement)
    [data.options.PsiGS0, data.options.PsiGS1, data.options.PsiGS2] = getDisplacements(data.options);
end;
clear lbpconfig imSize Images nImages gt gt_fname images_fname lbpconfig_fname options S

% create mapTable for extracting components of W
data.mapTable = f_createMapTable(data.options, data.tmpData);
toc
fprintf('Done.\n');

load('./MAT/kappasVAL.mat');

%% Validation

% load lambda values
load('./MAT/lambdas.mat');

%% start matlabpool
%try
%    matlabpool;
%catch ME1
%    idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Already connected to matlab pool. Skipping...\n');
%    end;
%end

fprintf('Running validation...\n');

Rval = zeros(1, length(lambdas));
% for lambda = lambdas
for j = 1 : length(lambdas)
    lambda = lambdas(j);
    load(['./results/W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '_lambda_' num2str(lambda) '.mat']);
    L = zeros(1, data.nImages);
    fprintf('VAL: Calculating Rval for lambda %f.\n', lambda);
    for par_i = 1 : data.nImages
        GT = data.Y{par_i};
        [Test TestSparse] = getPsiMat(data, par_i);
        P = argmax_mex(data.options, W, data.mapTable, TestSparse);        
        L(par_i) = 100*kappa(par_i) * 1/data.options.M * sum( sqrt( sum( (GT - P ).^2 ) ) );
    end;

    Rval(j) = 1/data.nImages * sum(L);
end;

%try
%    matlabpool close;
%catch ME2
%    idSegLast = regexp(ME2.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Cannot close matlabpool. Skipping...\n');
%    end;
%end;

fprintf('VAL: lambda = %s\n', num2str(lambdas, '%4.5f '));
fprintf('VAL: R      = %s\n', num2str(Rval, ' %4.5f '));

fprintf('Searching for minimum Rval...\n');

[mini idx] = min(Rval);
lambda = lambdas(idx);
fprintf('Minimum at lambda = %4.5f\n', lambda);

load(['./results/W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '_lambda_' num2str(lambda) '.mat']);
save(['./results/W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'], 'W');

fprintf('Optimal W saved.\n');

fprintf('Saving lambdas and Rvals... ');
save(['./results/Rvals_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'], 'lambdas', 'Rval');
fprintf('done.\n');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));
