%% test1.m - exp03
% This script evaluates vector W (see test_bmrm.m and testValidatoin.m)
% against test data set (TST).
% 
% 30-07-10 Michal Uricar
% 07-03-11 Michal Uricar, change to exp03

% clc; 
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Select functions

load('./MAT/functions_gdisp.mat');

%% Load TST data

% load config
load('./MAT/config.mat');
lbpconfig_fname = ['./MAT/lbpconfig_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
images_fname = ['./MAT/ImagesTst_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GTTst_nImages-' config.nImages '.mat'];

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

fprintf('Loading TST data...\n');
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

w_fname = ['./results/W_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
load(w_fname);

load('./MAT/kappas.mat');

%% Test

% % start matlabpool
%try
%    matlabpool;
%catch ME1
%    idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Already connected to matlab pool. Skipping...\n');
%    end;
%end

fprintf('Running test...\n');

err_nose = inf(1, data.nImages);
err_canthus_rl = inf(1, data.nImages);
err_canthus_lr = inf(1, data.nImages);
err_mouth_corner_r = inf(1, data.nImages);
err_mouth_corner_l = inf(1, data.nImages);
err_canthus_rr = inf(1, data.nImages);
err_canthus_ll = inf(1, data.nImages);
err_mean = inf(1, data.nImages);
err_maxdist = inf(1, data.nImages);
L = zeros(1, data.nImages);

tic

for i = 1 : data.nImages
    [Test TestSparse] = getPsiMat(data, i);
    GT = data.Y{i};
    P = argmax_mex(data.options, W, data.mapTable, TestSparse);
    
    errors = sqrt(sum((P - GT).^2)) * kappa(i) * 100;
    
    err_nose(i) = errors(1);
    err_canthus_rl(i) = errors(2);
    err_canthus_lr(i) = errors(3);
    err_mouth_corner_r(i) = errors(4);
    err_mouth_corner_l(i) = errors(5);
    err_canthus_rr(i) = errors(6);
    err_canthus_ll(i) = errors(7);
    err_mean(i) = mean(errors);
    err_maxdist(i) = max(errors);

    L(i) = 1/data.options.M * sum( errors );
end;
Rtst = 1/data.nImages * sum(L);
toc

%try
%    matlabpool close;
%catch ME2
%    idSegLast = regexp(ME2.identifier, '(?<=:)\w+$', 'match');
%    if (strcmp(idSegLast, 'OpenConnection'))
%        fprintf('Cannot close matlabpool. Skipping...\n');
%    end;
%end;

%% Save errors

fprintf('Saving errors...\n');
options = data.options;
save('./results/errors.mat', 'options', 'L', 'err_nose', 'err_canthus_rl', ...
            'err_canthus_lr', 'err_mouth_corner_r', 'err_mouth_corner_l', ...
            'err_canthus_rr', 'err_canthus_ll', 'err_mean', 'err_maxdist', 'Rtst');
fprintf('Done.\n');

   
fprintf('Rtst = %.4f\n', Rtst);

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));
