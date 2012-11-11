%% pruneDB.m
% This script prunes database to get rid of uncomplete or corrupted images
% 
% 15-07-10 Michal Uricar

% clc;
clearvars; close all;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Load old database

addpath('./Functions');

load('./MAT/eyefd_on_lfw.mat');

N = numel(image);

%% Prune database

% check if bad indices were already found and saved
if (~exist('./MAT/bad_indices.mat', 'file'))

    fprintf('Looking for bad indices...\n');
    image_path = './LFW/';
    indices = [];   

    for i = 1 : N
        itmp = image{i};

        filename = [image_path itmp.name];
        I = im2double(imread(filename));

        fname = strtok(itmp.name, '.');
        xmlfilename = [image_path 'mgt_v2/' fname '.xml'];
        T = face_XML_read(xmlfilename);

        fprintf('%.2f%% processing image %s...\n', i*100/N, itmp.name);

        if (~isfield(T, 'face'))
            fprintf('Incomplete annotation in file %s. Skipping...\n', filename);
            indices = [indices i];
            continue;
        end;

        cnt_bbox = size(itmp.bbox, 1);
        cnt_face = numel(T.face);

        if (cnt_bbox ~= 1)
            fprintf('Too much faces detected in file %s. Skipping...\n', filename);
            indices = [indices i];
            continue;
        end;

        if (cnt_face ~= 1)
            fprintf('Too much faces annotated in file %s. Skipping...\n', filename);
            indices = [indices i];
            continue;
        end;

        %if (~(isfield(T.face, 'eyel') && isfield(T.face, 'eyer') && ...
        %        isfield(T.face, 'mouth') && isfield(T.face, 'nose')))
        %    fprintf('Insufficient annotation in file %s. Skipping...\n', filename);
        %    indices = [indices i];
        %    continue;
        %end;
        if (~(isfield(T.face, 'canthus_rr') && isfield(T.face, 'canthus_rl') && ...
                isfield(T.face, 'canthus_lr') && isfield(T.face, 'canthus_ll') && ...
                isfield(T.face, 'mouth_corner_r') && isfield(T.face, 'mouth_corner_l') && ...
                isfield(T.face, 'nose') ))
            fprintf('Insufficient annotation in file %s. Skipping...\n', filename);
            indices = [indices i];
            continue;
        end;

    end;
    
    %% Save bad indices of original database

    fprintf('Saving bad indices...\n');
    save('./MAT/bad_indices.mat', 'indices');    
else
    fprintf('Loading bad indices...\n');
    load('./MAT/bad_indices.mat');
end;

%% Save pruned database

idx = 1:N;
idx(indices) = [];

image = image(idx);

fprintf('Saving pruned database...\n');
save('./MAT/eyefd_on_lfw_pruned.mat', 'image');

%% Timestamp

fprintf(1,'Finnished on %s\n\n', datestr(now));