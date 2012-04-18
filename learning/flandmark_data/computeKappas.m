%% computeKappas.m
% This scripts computes normalization coefficient for each image in database.
%
% 02-08-11 Michal Uricar

clc;
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% load GT - TRN

gt_fname = './MAT/GT_nImages-all.mat';
load(gt_fname);

%% compute kappas

kappa = inf(1, numel(gt));

for i = 1 : numel(gt)
    GT = gt{i};
    C = ( (GT(:, 2)+GT(:, 6))/2 + (GT(:, 3)+GT(:, 7))/2 )/2;
    kappa(i) = 1/sqrt(sum(( C - (GT(:, 4)+GT(:, 5))/2 ).^2));
end;

%% save kappas

save('./MAT/kappas.mat', 'kappa');

%% load GT - VAL

gt_fname = './MAT/GTVal_nImages-all.mat';
load(gt_fname);

%% compute kappas

kappa = inf(1, numel(gt));

for i = 1 : numel(gt)
    GT = gt{i};
    C = ( (GT(:, 2)+GT(:, 6))/2 + (GT(:, 3)+GT(:, 7))/2 )/2;
    kappa(i) = 1/sqrt(sum(( C - (GT(:, 4)+GT(:, 5))/2 ).^2));
end;

%% save kappas

save('./MAT/kappasVAL.mat', 'kappa');

%% load GT - TST

gt_fname = './MAT/GTTst_nImages-all.mat';
load(gt_fname);

%% compute kappas

kappa = inf(1, numel(gt));

for i = 1 : numel(gt)
    GT = gt{i};
    C = ( (GT(:, 2)+GT(:, 6))/2 + (GT(:, 3)+GT(:, 7))/2 )/2;
    kappa(i) = 1/sqrt(sum(( C - (GT(:, 4)+GT(:, 5))/2 ).^2));
end;

%% save kappas

save('./MAT/kappasTST.mat', 'kappa');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));

