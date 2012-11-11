%% dbDivision.m
% This script randomly divides database into 3 parts - training (TRN),
% validating (VAL) and testing (TST).
% 
% 16-07-10 Michal Uricar
% 21-03-12 Michal Uricar, LFW annotation in one file

clearvars; close all;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Load pruned database

addpath('./Functions');

% load('./MAT/eyefd_on_lfw_pruned.mat');
load('./MAT/db_good.mat');

tmp_as = annotation_struct;

N = annotation_struct.N;

per60 = round(60/100*N);
per20 = round(20/100*N);

tic
indices_perm = randperm(N);
toc

% 60% of images will form training set
training = indices_perm(1:per60); indices = training;
indices_perm(1:per60) = []; 
% save training set
%image = tmp_image(training);
annotation_struct.bbox = tmp_as.bbox(training, :);
annotation_struct.names = tmp_as.names(training);
annotation_struct.eye_r = tmp_as.eye_r(training, :);
annotation_struct.eye_l = tmp_as.eye_l(training, :);
annotation_struct.canthus_rr = tmp_as.canthus_rr(training, :);
annotation_struct.canthus_rl = tmp_as.canthus_rl(training, :);
annotation_struct.canthus_lr = tmp_as.canthus_lr(training, :);
annotation_struct.canthus_ll = tmp_as.canthus_ll(training, :);
annotation_struct.mouth = tmp_as.mouth(training, :);
annotation_struct.mouth_corner_r = tmp_as.mouth_corner_r(training, :);
annotation_struct.mouth_corner_l = tmp_as.mouth_corner_l(training, :);
annotation_struct.nose = tmp_as.nose(training, :);
annotation_struct.N = numel(annotation_struct.names);

% save('./MAT/training.mat', 'image');
save('./MAT/TRN.mat', 'annotation_struct', 'indices');

% 20% of images will form validating set
validating = indices_perm(1:per20); indices = validating;
indices_perm(1:per20) = [];
% save validating set
% image = tmp_image(validating);
annotation_struct.bbox = tmp_as.bbox(validating, :);
annotation_struct.names = tmp_as.names(validating);
annotation_struct.eye_r = tmp_as.eye_r(validating, :);
annotation_struct.eye_l = tmp_as.eye_l(validating, :);
annotation_struct.canthus_rr = tmp_as.canthus_rr(validating, :);
annotation_struct.canthus_rl = tmp_as.canthus_rl(validating, :);
annotation_struct.canthus_lr = tmp_as.canthus_lr(validating, :);
annotation_struct.canthus_ll = tmp_as.canthus_ll(validating, :);
annotation_struct.mouth = tmp_as.mouth(validating, :);
annotation_struct.mouth_corner_r = tmp_as.mouth_corner_r(validating, :);
annotation_struct.mouth_corner_l = tmp_as.mouth_corner_l(validating, :);
annotation_struct.nose = tmp_as.nose(validating, :);
annotation_struct.N = numel(annotation_struct.names);
% save('./MAT/validating.mat', 'image');
save('./MAT/VAL.mat', 'annotation_struct', 'indices');

% 20% (the rest) of images will form testing set
testing = indices_perm; indices = testing;
% save testing set
% image = tmp_image(testing);
annotation_struct.bbox = tmp_as.bbox(testing, :);
annotation_struct.names = tmp_as.names(testing);
annotation_struct.eye_r = tmp_as.eye_r(testing, :);
annotation_struct.eye_l = tmp_as.eye_l(testing, :);
annotation_struct.canthus_rr = tmp_as.canthus_rr(testing, :);
annotation_struct.canthus_rl = tmp_as.canthus_rl(testing, :);
annotation_struct.canthus_lr = tmp_as.canthus_lr(testing, :);
annotation_struct.canthus_ll = tmp_as.canthus_ll(testing, :);
annotation_struct.mouth = tmp_as.mouth(testing, :);
annotation_struct.mouth_corner_r = tmp_as.mouth_corner_r(testing, :);
annotation_struct.mouth_corner_l = tmp_as.mouth_corner_l(testing, :);
annotation_struct.nose = tmp_as.nose(testing, :);
annotation_struct.N = numel(annotation_struct.names);
% save('./MAT/testing.mat', 'image');
save('./MAT/TST.mat', 'annotation_struct', 'indices');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));