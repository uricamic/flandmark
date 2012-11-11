%% runstats.m
% This script collects statistics about positions of each component,
% creates matrix S 4xM [xmin_0, ..., xmin_{M-1}; ymin_0, ..., ymin_{M-1};
% xmin_0, ..., xmax_{M-1}; ymax_0, ..., ymax_{M-1}] each column defines
% bounding box for one component. 
% 
% Parameters in structure options are now set to cover approx. 97.5% of all
% tested images (about 12200 images, original database was pruned - file
% ./MAT/eye_on_lfw_pruned.mat - so now it does not contain images with
% multiple faces annotated or detected. Also images with incomplete
% annotations were discarded.) 
% 
% 15-07-10 Michal Uricar
% 26-02-11 Michal Uricar

% clc;
clearvars; close all;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Load data (images with detected frames)

addpath('./Functions');
addpath('../../matlab_toolbox/mex/');

load('./MAT/paths.mat');
load('./MAT/options.mat');

% load('./MAT/eyefd_on_lfw_pruned.mat');
annotation_struct = read_lfw_annotation_file(['../' DB_ANN]);

%% Print parameters

fprintf('Parameters: \nbw = [%d; %d]\n', options.bw);
fprintf('bw_margin = [%d; %d]\n', options.bw_margin);
fprintf('components= \n'); disp(options.components);
fprintf('M = %d\n', options.M);
fprintf('image_path = ''%s''\n', options.image_path);
fprintf('....................................\n');
fprintf('\nStarting test...\n');

%% run test

tic
[cnt, S, bad_idx] = paramStats(annotation_struct, options);
% [cnt, S, bad_idx] = paramStats(annotation_struct, options, 1);
toc

percent = cnt/annotation_struct.N * 100;
fprintf('Passed %d images. That is %f%%.\n', cnt, percent);

%% create new database - only images that passed through test

idx = 1:annotation_struct.N;
idx(bad_idx) = [];
% image = image(idx);
annotation_struct.bbox = annotation_struct.bbox(idx, :);
annotation_struct.names = annotation_struct.names(idx);
annotation_struct.eye_r = annotation_struct.eye_r(idx, :);
annotation_struct.eye_l = annotation_struct.eye_l(idx, :);
annotation_struct.canthus_rr = annotation_struct.canthus_rr(idx, :);
annotation_struct.canthus_rl = annotation_struct.canthus_rl(idx, :);
annotation_struct.canthus_lr = annotation_struct.canthus_lr(idx, :);
annotation_struct.canthus_ll = annotation_struct.canthus_ll(idx, :);
annotation_struct.mouth = annotation_struct.mouth(idx, :);
annotation_struct.mouth_corner_r = annotation_struct.mouth_corner_r(idx, :);
annotation_struct.mouth_corner_l = annotation_struct.mouth_corner_l(idx, :);
annotation_struct.nose = annotation_struct.nose(idx, :);
annotation_struct.N = numel(annotation_struct.names);

fprintf('Creating new database (db_good.mat)...\n');
save('./MAT/db_good.mat', 'annotation_struct');

%% Save data

fprintf('Saving statistics...\n');
save('./results/statistics_orig.mat', 'cnt', 'S', 'percent', 'options');

 %% shrink S (S boundary must represent positions of components center)

load('./results/statistics_orig.mat', 'S');
 
for i = 1 : options.M
    S(1, i) = S(1, i) + options.components(1, i)/2;
    S(2, i) = S(2, i) + options.components(2, i)/2;
    S(3, i) = S(3, i) - options.components(1, i)/2;
    S(4, i) = S(4, i) - options.components(2, i)/2;
end;

save('./results/statistics_orig.mat', 'cnt', 'S', 'percent', 'options');

%% Show computed S

[Iframe, Annotation] = getImageFrame(options, 1, annotation_struct);
Points = Annotation.P;
Points = prepareS0gt(Annotation.P, options);    % transform nose to the center of face
Points = Points(:, options.comselect);          % extract relevant points only (name list in options.compnames)
Points(:, options.M) = Annotation.P(:, 10);             % copy back original nose position

colors = colormap(hsv(options.M)); close gcf;
% names = fieldnames(Annotation.face);
names = options.compnames;

scrsz = get(0,'ScreenSize');

figure;
imshow(Iframe, [], 'Border', 'tight'); hold on;
set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
% set(gcf, 'Position', [100, 100, 800, 600]);
for k = 1 : options.M
%     plot(Points(1, k), Points(2, k), 'r.', 'LineWidth', 2, 'MarkerSize', 10);
    plot(Points(1, k), Points(2, k), '.', 'color', colors(k, :), 'LineWidth', 2, 'MarkerSize', 10);
    bb = makeAABB(S(:, k));
    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', colors(k, :));
    bb = makeAABB([Points(1, k) - options.components(1, k)/2 Points(2, k) - options.components(2, k)/2 ...
                   Points(1, k) + options.components(1, k)/2 Points(2, k) + options.components(2, k)/2 ]);
%     text(Points(1, k)-1, Points(2, k)+1, names(k+3), 'color', colors(k, :));
    text(Points(1, k)-1, Points(2, k)+1, names(k), 'color', colors(k, :));
    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', 'y');
end;
hold off;

if (~exist('./img/', 'dir'))
    mkdir('./img');
end;

saveas(gcf, './img/S.png');

%% Show S

S = ceil(S);

colors = colormap(hsv(options.M)); close gcf;
% names = fieldnames(Annotation.face);
names = options.compnames;
i = 1;

figure;
subplot(1, 2, 1);
imshow(Iframe, []); hold on;
for j = 1 : options.M;
    bb = makeAABB(S(:, j));
    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', colors(j, :));
    text(Points(1, j), Points(2, j), names(j), 'color', colors(j, :));
end;
line([20 20], [0 41], 'color', 'm', 'LineWidth', 2);
hold off;
saveas(gcf, './img/S.png');

%% make each frame of component in S symetric (this works only on lfw database with 4 components!)

S_new = S;

% canthus - inner
S_new(2, 2) = min(S(2, 2), S(2, 3)); S_new(2, 3) = min(S(2, 2), S(2, 3));   % ymin
S_new(4, 2) = max(S(4, 2), S(4, 3)); S_new(4, 3) = max(S(4, 2), S(4, 3));   % ymax
S_new(3, 2) = options.bw(1)/2 + max(S(3, 2) - options.bw(1)/2, options.bw(1)/2 - S(1, 3));     % x_max
S_new(1, 3) = options.bw(1)/2 - max(S(3, 2) - options.bw(1)/2, options.bw(1)/2 - S(1, 3));
S_new(1, 2) = min(S(1, 2), options.bw(1) - S(3, 3));
S_new(3, 3) = options.bw(1) - min(S(1, 2), options.bw(1) - S(3, 3));
% canthus - outer
S_new(2, 6) = min(S(2, 6), S(2, 7)); S_new(2, 7) = min(S(2, 6), S(2, 7));   % ymin
S_new(4, 6) = max(S(4, 6), S(4, 7)); S_new(4, 7) = max(S(4, 6), S(4, 7));   % ymax
S_new(3, 6) = options.bw(1)/2 + max(S(3, 6) - options.bw(1)/2, options.bw(1)/2 - S(1, 7));     % x_max
S_new(1, 7) = options.bw(1)/2 - max(S(3, 6) - options.bw(1)/2, options.bw(1)/2 - S(1, 7));
S_new(1, 6) = min(S(1, 6), options.bw(1) - S(3, 7));
S_new(3, 7) = options.bw(1) - min(S(1, 6), options.bw(1) - S(3, 7));
% mouth corners
S_new(2, 4) = min(S(2, 4), S(2, 5)); S_new(2, 5) = min(S(2, 4), S(2, 5));   % ymin
S_new(4, 4) = max(S(4, 4), S(4, 5)); S_new(4, 5) = max(S(4, 4), S(4, 5));   % ymax
S_new(3, 4) = options.bw(1)/2 + max(S(3, 4) - options.bw(1)/2, options.bw(1)/2 - S(1, 5));     % x_max
S_new(1, 5) = options.bw(1)/2 - max(S(3, 4) - options.bw(1)/2, options.bw(1)/2 - S(1, 5));
S_new(1, 4) = min(S(1, 4), options.bw(1) - S(3, 5));
S_new(3, 5) = options.bw(1) - min(S(1, 4), options.bw(1) - S(3, 5));
% center of face
S_new(3, 1) = options.bw(1)/2 + max(S(3, 1) - options.bw(1)/2, options.bw(1)/2 - S(1, 1));     % x_max
S_new(1, 1) = options.bw(1)/2 - max(S(3, 1) - options.bw(1)/2, options.bw(1)/2 - S(1, 1));
% nose
S_new(3, 8) = options.bw(1)/2 + max(S(3, 8) - options.bw(1)/2, options.bw(1)/2 - S(1, 8));     % x_max
S_new(1, 8) = options.bw(1)/2 - max(S(3, 8) - options.bw(1)/2, options.bw(1)/2 - S(1, 8));

subplot(1, 2, 2);
imshow(Iframe); hold on;
for j = 1 : options.M;
    bb = makeAABB(S_new(:, j));
    line([bb(1,:) bb(1,1)], [bb(2,:) bb(2,1)], 'color', colors(j, :));
    text(Points(1, j), Points(2, j), names(j), 'color', colors(j, :));
end;
line([20 20], [0 41], 'color', 'm', 'LineWidth', 2);
hold off;
saveas(gcf, './img/S_symetric.png');

%% Save symetric S for 4 components
S2 = S_new;
% save('./results/symetric_S_orig.mat', 'S');

% set new options.S
options.S = S2;

fname = './results/options.mat';
fprintf('Saving new options struct to file %s...', fname);
save(fname, 'options');
fprintf(' done.\n');

%% Check for TRN, VAL and TST annotation files

% if (exist('./MAT/LFW_TRN_annotation.txt', 'file') && exist('./MAT/LFW_VAL_annotation.txt', 'file') && exist('./MAT/LFW_TST_annotation.txt', 'file'))
%     annotation_struct = read_lfw_annotation_file('./MAT/LFW_TRN_annotation.txt');
%     save('./MAT/TRN.mat', 'annotation_struct');
%     annotation_struct = read_lfw_annotation_file('./MAT/LFW_VAL_annotation.txt');
%     save('./MAT/VAL.mat', 'annotation_struct');
%     annotation_struct = read_lfw_annotation_file('./MAT/LFW_TST_annotation.txt');
%     save('./MAT/TST.mat', 'annotation_struct');
% end;
if (~isempty(TRN_EXAMPLES)  && ~isempty(VAL_EXAMPLES) && ~isempty(TST_EXAMPLES))
    if (exist(['../' TRN_EXAMPLES], 'file') && exist(['../' VAL_EXAMPLES], 'file') && exist(['../' TST_EXAMPLES], 'file'))
        annotation_struct = read_lfw_annotation_file(['../' TRN_EXAMPLES]);
        save('./MAT/TRN.mat', 'annotation_struct');
        annotation_struct = read_lfw_annotation_file(['../' VAL_EXAMPLES]);
        save('./MAT/VAL.mat', 'annotation_struct');
        annotation_struct = read_lfw_annotation_file(['../' TST_EXAMPLES]);
        save('./MAT/TST.mat', 'annotation_struct');
    end;
end;

%% Timestamp

fprintf(1,'\nFinished on %s\n\n', datestr(now));