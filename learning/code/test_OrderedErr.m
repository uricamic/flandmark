%% test_OrderedErr.m - exp03
% This script automatically generates HTML page with images sorted
% according to error.
% 
% 28-10-10 Michal Uricar
% 01-03-11 Michal Uricar

clc;
clearvars; close all;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');
addpath('../../matlab_toolbox/mex/');

%% Options

% show images
opt.verbose = true;
% save images
% opt.save = true;
opt.save = false;

%% Select functions

load('./MAT/functions_gdisp.mat');
gfunc = 'gdisp';

load('./results/W_nImages-all_hop-4.mat');
model.W = W;
model.fgetPsiMat = @getPsiMat;
model.fgetQG = @getQG;
model.fargmax = @argmax2;

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

data.lbp = lbpconfig; data.imSize = imSize; data.Images = Images; data.nImages = nImages;
options.S = S; data.Y = gt; data.options = options;
data.tmpData = getPsiMat(data, 1);
if (displacement)
 [data.options.PsiGS0, data.options.PsiGS1, data.options.PsiGS2] = getDisplacements(data.options);
end;
data.mapTable = f_createMapTable(data.options, data.tmpData);

model.data = data;

clearvars -except model gfunc;

save('detector_model.mat', 'model');

%% Load data

% select dataset ['TST','TRN','VAL']
dataset = 'TST';

switch dataset
    case 'TRN'
        images_fname = './MAT/Images_nImages-all.mat';
        gt_fname = './MAT/GT_nImages-all.mat';
    case 'VAL'
        images_fname = './MAT/ImagesVal_nImages-all.mat';
        gt_fname = './MAT/GTVal_nImages-all.mat';
    case 'TST'
        images_fname = './MAT/ImagesTst_nImages-all.mat';
        gt_fname = './MAT/GTTst_nImages-all.mat';
    otherwise 
        fprintf('Fatal Error: dataset must be one of [''TST'', ''TRN'', ''VAL'']\n');
        exit;
end;

load(images_fname);
load(gt_fname);
load(['./MAT/' dataset '.mat']);

load('./MAT/kappas.mat');

%% Generate HTML page with error-ordered TST images

if (~exist('./html', 'dir'))
    mkdir('./html');
end;
if (~exist(['./html/img/' gfunc '_' dataset], 'dir'))
    mkdir(['./html/img/' gfunc '_' dataset]);
end;

html_file = fopen(['./html/' gfunc '_' dataset '.php'], 'w');

%% Compute detection

if (~exist(['./html/errors' dataset '_' gfunc '.mat'], 'file'))
    nImages = numel(image);

    err_nose = inf(1, nImages);
    err_canthus_rl = inf(1, nImages);
    err_canthus_lr = inf(1, nImages);
    err_mouth_corner_r = inf(1, nImages);
    err_mouth_corner_l = inf(1, nImages);
    err_canthus_rr = inf(1, nImages);
    err_canthus_ll = inf(1, nImages);
    err_mean = inf(1, nImages);
    err_maxdist = inf(1, nImages);
    L = zeros(1, nImages);

    S_max = cell(1, nImages);

    comps = ['S0'; 'S1'; 'S2'; 'S3'; 'S4'; 'S5'; 'S6'];

    for i = 1 : nImages
        fname = [image{i}.name(1:end-4) '.png'];

        %[S_max{i}, Q{i}, G{i}] = detector(Images(:, i), model);
	[S_max{i}] = detector(Images(:, i), model);

        GT = gt{i};
        I = reshape(Images(:, i), imSize);

        Smax = S_max{i};

        r = 251; c = 251;
        fig = figure('Visible', 'off');
        imshow(I, [], 'Border', 'tight'); hold on;
        plot(GT(1, :), GT(2, :), 'gs', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'g');
        text(GT(1, :), GT(2, :)+1, comps, 'color', 'g', 'FontSize', 12);
        plot(Smax(1, :), Smax(2, :), 'rs', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        text(Smax(1, :), Smax(2, :)+1, comps, 'color', 'r', 'FontSize', 12);
        hold off;

        % Sets position and size of figure on the screen
        set(fig, 'Units', 'pixels', 'position', [250 250 c r] ); 
        % Sets axes to fill the figure space
        set(gca, 'Units', 'pixels', 'position', [0 0 c+1 r+1 ]);
        % Sets print properties; Looks like 1 pixel = (3/4)th of a point
        set(fig, 'paperunits', 'points', 'papersize', [fix((c-1)*(3/4))+1 fix((r-1)*(3/4))+1]);
        set(fig, 'paperunits', 'normalized', 'paperposition', [0 0 1 1]);
        print( fig, sprintf('-r%d', ceil(72*(4/3))), '-dpng', ['./html/img/' gfunc '_' dataset '/' fname] );

        close all;

        %errors = sqrt(sum((S_max{i} - GT).^2));
	errors = sqrt(sum((S_max{i} - GT).^2)) * kappa(i) * 100;

        % order of components S0, S1, S2, S3 - 'nose', 'eyel', 'mouth', 'eyer'
        err_nose(i) = errors(1);
        err_canthus_rl(i) = errors(2);
        err_canthus_lr(i) = errors(3);
        err_mouth_corner_r(i) = errors(4);
        err_mouth_corner_l(i) = errors(5);
        err_canthus_rr(i) = errors(6);
        err_canthus_ll(i) = errors(7);
        err_mean(i) = mean(errors);
        err_maxdist(i) = max(errors);

        L(i) = 1/model.data.options.M * sum( errors );

        fprintf('Image %d done ... %f%%\n', i, i/nImages*100);
    end;
    R = 1/nImages * sum(L);

    %save(['./html/errors' dataset '_' gfunc '.mat'], 'model', 'L', 'err_eyel', 'err_eyer', 'err_mouth', 'err_nose',...
    %    'err_mean', 'err_maxdist', 'R', 'S_max', 'Q', 'G', 'comps');
    save(['./html/errors' dataset '_' gfunc '.mat'], 'L', 'err_nose', 'err_canthus_rl', ...
            'err_canthus_lr', 'err_mouth_corner_r', 'err_mouth_corner_l', ...
            'err_canthus_rr', 'err_canthus_ll', 'err_mean', 'err_maxdist', 'R', 'comps');

else
    load(['./html/errors' dataset '_' gfunc '.mat']);
end

%% generate HTML file with ordered errors

cnt_errpx = [];
pos = 1;
from = 0;

for errpx = 0 : model.data.options.bw(1)    
    
    [err_mean_srt, srt_idx] = sort(err_mean);
    meanerr = round(err_mean);
    A = find(meanerr == errpx);
    to = from+length(A);
    if (isempty(A))
        best_inl = [];
    else
        best_inl = srt_idx(from+1:to);
    end;
    from = to;
        
    cnt_errpx(1, pos) = errpx;
    cnt_errpx(2, pos) = length(best_inl);
    pos = pos + 1;
    
    if (length(best_inl) > 0)
        fprintf(html_file, ...
            [...
            '<h3><a name="errpx' num2str(errpx) '">Errpx ' num2str(errpx) '</a></h3>\n',...
            '\n<table cellspacing="0">\n',...
            '<caption>Measured erros (' gfunc ')</caption>\n',...
            '<tr>\n',...
            '<th scope="col">Image</th>\n',...
            '<th scope="col">MeanErr</th>\n',...
            '<th scope="col">MaxErr</th>\n',...
            ]);

        for j = 0 : model.data.options.M-1
            fprintf(html_file, ['<th scope="col">q(I, S_' num2str(j) ')</th>\n']);
        end;
        fprintf(html_file, '<th scope="col">sum_{i=1}^{M-1}{q(I, S_i)}</th>\n');
        for j = 1 : model.data.options.M-1
            fprintf(html_file, ['<th scope="col">g(S_0, S_' num2str(j) ')</th>\n']);
        end;
        fprintf(html_file, '<th scope="col">sum_{i=1}^{M-1}{g(S_0, S_i)}</th>\n');

        % ordered error
        for i = 1 : length(best_inl)
            idx = best_inl(i);
            fname = [image{idx}.name(1:end-4) '.png'];
            Qi = Q{idx}; Gi = G{idx};

            fprintf(html_file, ['<tr>\n<td><img src="' './img/' gfunc '_' dataset '/' fname '" /></td>\n']);
            fprintf(html_file, ['<td>' num2str(err_mean(idx)) '</td>\n<td>' num2str(err_maxdist(idx)) '</td>\n']);
            
            for j = 1 : model.data.options.M
                fprintf(html_file, ['<td>' num2str(Qi(j)) '</td>\n']);
            end;
            fprintf(html_file, ['<td class="sum">' num2str(sum(Qi)) '</td>\n']);
            for j = 1 : model.data.options.M - 1
                fprintf(html_file, ['<td>' num2str(Gi(j)) '</td>\n']);
            end;
            fprintf(html_file, ['<td class="sum">' num2str(sum(Gi)) '</td>\n']);
            fprintf(html_file, '</tr>\n');

            fprintf('Errpx %d, Done ... %f%%\n', errpx, i/length(best_inl)*100);
        end;

        fprintf(html_file, '</table>\n');
    end;
end;

cnt_errpx(:, cnt_errpx(2, :)==0) = [];

fprintf('Errpx: \t');
for i = 1 : size(cnt_errpx, 2)
    fprintf('%d \t', cnt_errpx(1, i));
end;
fprintf('\n');
fprintf('Cnt: \t');
for i = 1 : size(cnt_errpx, 2)
    fprintf('%d \t', cnt_errpx(2, i));
end;
fprintf('\n');

save('./html/cnt_errpx.mat', 'cnt_errpx');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));