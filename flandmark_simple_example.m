%% flandmark_simple_example.m
% Test flandmark detector
% 
% 27-03-12 Michal Uricar

clc; 
clearvars; close all;

%% Add path

addpath('./learning/flandmark_data/Functions/');
addpath('./learning/flandmark_data/libocas/mex/');

DIR = './data/Images/';
IMGS = dir([DIR '*.jpg']);

% pick random image
idx = randi(numel(IMGS));
filename = IMGS(idx).name;

tokens = strtokenize(filename, '.');
fname = tokens(1:end-1);
fname = fname{1};

% load image and detected face bbox
I = imread([DIR filename]);
bbox = dlmread([DIR fname '.det']);

% image output
figure(1);
imshow(I, [], 'Border', 'tight'); hold on;
plotbox(bbox);

for i = 1 : size(bbox, 1)
    tic
    P = flandmark_detector(rgb2gray(I), int32(bbox(i, :)),  './data/flandmark_model.dat');
    t1 = toc;
    fprintf('MEX:    Elapsed time %f ms\n', t1*1000);

    % show landmarks 
    comps = ['S0'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6'; 'S7'];
    plot(P(1, 1), P(2, 1), 'ys', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'y');
    text(P(1, 1)+1, P(2, 1)+1, comps(1,:), 'color', 'y', 'FontSize', 12);
    plot(P(1, 2:end), P(2, 2:end), 'rs', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    text(P(1, 2:end)+1, P(2, 2:end)+1, comps(2:end,:), 'color', 'r', 'FontSize', 12);
end;
    
clear mex;