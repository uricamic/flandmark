%% detectorModelBuild.m
% This script creates mat file with detector's model structure (for use
% with the flandmark_detector.
% 
% 29-07-11 Michal Uricar

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Create mat file

load('./MAT/functions_gdisp.mat');
gfunc = 'gdisp';

load('./results/W_nImages-all_hop-4.mat');
model.W = W;
model.fgetPsiMat = @getPsiMat;
model.fgetQG = @getQG;
model.fargmax = @argmax2;
model.argmax_mex = @argmax_mex;

load('./MAT/config.mat');

lbpconfig_fname = './MAT/lbpconfig_nImages-all_hop-4.mat';
images_fname = './MAT/ImagesTst_nImages-all.mat';
gt_fname = './MAT/GTTst_nImages-all.mat';

load(lbpconfig_fname);
load(images_fname);
load(gt_fname);

load('./results/options.mat');

data.lbp = lbpconfig; 
data.imSize = imSize; 
data.Images = Images(:, 1); 
data.nImages = size(Images, 2);
data.Y = gt; data.options = options;
data.tmpData = getPsiMat(data, 1);
[data.options.PsiGS0, data.options.PsiGS1, data.options.PsiGS2] = getDisplacements(data.options);
data.mapTable = f_createMapTable(data.options, data.tmpData);
model.data = data;
clearvars -except model gfunc;

save('model.mat', 'model');

%% Generate binary .dat file
% 
% addpath('./flandmark_mex/');

flandmark_write_model(model, './flandmark_model.dat');

clear mex;

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));
