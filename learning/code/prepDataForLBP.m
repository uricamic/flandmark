%% prepDataForLBP.m
% Loads prepared datasets (TRN, VAL, TST) and prepares data for LBP
% features computation.
% 
% 12-07-11 Michal Uricar

clc; 
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% add paths

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Load data

overwrite = false;
% overwrite = true;   % re-generate images and gt mat files

% load options
load('./results/options.mat');
S = options.S;
% verbose output
options.verbose = true;

M = size(options.components, 2);

% load config
load('./MAT/config.mat');

%% Load all images & save annotations

% TRN
images_fname = ['./MAT/Images_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GT_nImages-' config.nImages '.mat'];
if (~exist(images_fname, 'file') || overwrite)
    load('./MAT/TRN.mat');
    bad_idx_TRN = prepLBPconf(config, options, 'TRN', annotation_struct, images_fname, gt_fname);
else
    fprintf('TRN: Array of images already created. Loading %s...\n', images_fname);
    load(images_fname);
end;

% VAL
images_fname = ['./MAT/ImagesVal_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GTVal_nImages-' config.nImages '.mat'];
if (~exist(images_fname, 'file') || overwrite)
    load('./MAT/VAL.mat');
    bad_idx_VAL = prepLBPconf(config, options, 'VAL', annotation_struct, images_fname, gt_fname);
else
    fprintf('VAL: Array of images already created. Loading %s...\n', images_fname);
    load(images_fname);
end;

% TST
images_fname = ['./MAT/ImagesTst_nImages-' config.nImages '.mat'];
gt_fname = ['./MAT/GTTst_nImages-' config.nImages '.mat'];
if (~exist(images_fname, 'file') || overwrite)
    load('./MAT/TST.mat');
    bad_idx_TST = prepLBPconf(config, options, 'TST', annotation_struct, images_fname, gt_fname);
else
    fprintf('TST: Array of images already created. Loading %s...\n', images_fname);
    load(images_fname);
end;

if (exist('bad_idx_TRN', 'var') && exist('bad_idx_VAL', 'var') && exist('bad_idx_TST', 'var'))
    save('./MAT/bad_idx.mat', 'bad_idx_TRN', 'bad_idx_VAL', 'bad_idx_TST');
end;

%% Prepare data for LBP

fprintf('Preparing data for LBP...\n');

lbpconfig = cell(1, M);
    
% get psi matrix for each component
for j = 1 : M
    siz = [S(4, j) - S(2, j) + 1; S(3, j) - S(1, j) + 1];
    offset = S(1:2, j);
    coors = zeros(2, siz(1)*siz(2));
    for cor_i = 1 : siz(2)
        coors(1, ((cor_i-1)*siz(1)+1):cor_i*siz(1)) = (cor_i - 1 + offset(1));
        coors(2, ((cor_i-1)*siz(1)+1):cor_i*siz(1)) = (offset(2):(siz(1)-1+offset(2)));
    end;
    
    card_si = size(coors, 2);
    
    height = round(options.components(2, j));
    width = round(options.components(1, j));
    lbpconfig{j}.winSize = [height; width];
    lbpconfig{j}.hop = config.heightOfPyramid;
    lbpconfig{j}.wins = uint32([ones(1, card_si); coors; zeros(1, card_si)] - repmat([0; floor(width/2); floor(height/2); 0], 1, card_si));
end;

lbpconfig_fname = ['./MAT/lbpconfig_nImages-' config.nImages '_hop-' num2str(config.heightOfPyramid) '.mat'];
fprintf('Saving data %s ... ', lbpconfig_fname);
save(lbpconfig_fname, 'lbpconfig');
% save(lbpconfig_fname, 'Images', 'gt', 'lbpconfig', 'options', 'S');
fprintf('Done.\n');
