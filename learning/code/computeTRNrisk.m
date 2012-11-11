%% computeTRNrisk.m
% This script does validation of vectors W acquired from training data
% (see test_bmrm.m)
% 
% 27-04-11 Michal Uricar

clc; 
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Select functions

load('./MAT/functions_gdisp.mat');

%% 
% load config
load('./MAT/config.mat');
lbpconfig_fname = ['./MAT/exp03_lbpconfig_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
images_fname = ['./MAT/Images_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GT_nImages-' config.nImages '.mat'];

fprintf('Loading VAL data...\n');
load(lbpconfig_fname);
load(images_fname);
load(gt_fname);
load('./results/exp03_options.mat');
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
    data.options.PsiG = getDisplacements(data.options);
end;
clear lbpconfig imSize Images nImages gt gt_fname images_fname lbpconfig_fname options S

% create mapTable for extracting components of W
data.mapTable = f_createMapTable(data.options, data.tmpData);
toc
fprintf('Done.\n');

%% 

lambdas = [1e-2, 1e-1, 1];

% start matlabpool
try
    matlabpool;
catch ME1
    idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
    if (strcmp(idSegLast, 'OpenConnection'))
        fprintf('Already connected to matlab pool. Skipping...\n');
    end;
end

fprintf('Running validation...\n');

Rtrn = zeros(1, length(lambdas));
% for lambda = lambdas
for j = 1 : length(lambdas)
    lambda = lambdas(j);
    load(['./results/exp03_W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '_lambda_' num2str(lambda) '.mat']);
    L = zeros(1, data.nImages);
    fprintf('TRN: Calculating Rtrn for lambda %f.\n', lambda);
%     for i = 1 : data.nImages
    parfor par_i = 1 : data.nImages
        GT = data.Y{par_i};
        Test = getPsiMat(data, par_i);
        [ q, g ] = getQG(data.options, W, Test, data.mapTable);
        P = f_argmax(data.options, q, g);
        L(par_i) = 1/data.options.M * sum( sqrt( sum( (GT - P ).^2 ) ) );
    end;
    Rtrn(j) = 1/data.nImages * sum(L);
end;

try
    matlabpool close;
catch ME2
    idSegLast = regexp(ME2.identifier, '(?<=:)\w+$', 'match');
    if (strcmp(idSegLast, 'OpenConnection'))
        fprintf('Cannot close matlabpool. Skipping...\n');
    end;
end;

fprintf('TRN: lambda = %s\n', num2str(lambdas, '%4.5f '));
fprintf('TRN: R      = %s\n', num2str(Rtrn, ' %4.5f '));

fprintf('Saving lambdas and Rtrns... ');
save(['./results/exp03_Rtrns_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'], 'lambdas', 'Rtrn');
fprintf('done.\n');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));
