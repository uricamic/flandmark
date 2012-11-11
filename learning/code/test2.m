%% test2.m - exp03
% This script shows histograms for errors computed in test1.m
% 
% 02-08-10 Michal Uricar
% 07-03-11 Michal Uricar

% clc;
close all; clearvars;

%% Timestamp

fprintf(1,'Started on %s\n\n', datestr(now));

%% Add path

addpath('./Functions/');

%% Load errors

load('./results/errors.mat');

%% Options

% show images
% opt.verbose = true;
opt.verbose = false;
% save images
opt.save = true;
% opt.save = false;

%% Show histograms

if (~exist('./img/', 'dir'))
    mkdir('./img/');
end;

% get screen size
scrsz = get(0,'ScreenSize');

% histogram of mean errors
x = 0:options.bw(1); y = round(L);
hist_mean = hist(y, x);
% cumulative histogram of mean errors
n_elements = histc(y, x); c_elements = cumsum(n_elements);
if (opt.verbose)
    figure; hist(y, x);
    set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    set(gca, 'XTick', x); set(gca, 'xlim', [-0.5 options.bw(1)+0.5]);
    xlabel('Distance in pixels'); ylabel('Count of occurrences');
    set(gca, 'YTick', 0:50:2*c_elements(end));
    title('Histogram of mean errors');
    if (opt.save)
        saveas(gcf, './img/histogram_mean.png');
        saveas(gcf, './img/histogram_mean.fig');
    end;
    
    figure; bar(x, c_elements / max(c_elements));
    xlabel('RMS [px]'); ylabel('Count of occurrences [%]');
    set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    set(gca, 'XTick', x); set(gca, 'xlim', [-0.5 options.bw(1)+0.5]);
    set(gca, 'YTick', 0:0.05:1);
    title('Cumulative historgam of mean errors');
    if (opt.save)
        saveas(gcf, './img/cumhist_mean.png');
        saveas(gcf, './img/cumhist_mean.fig');
    end;
end;


% histogram of max errors
x = 0:options.bw(1); y = round(err_maxdist);
hist_max = hist(y, x);
% cumulative histogram of max errors
n_elements = histc(y, x); c_elements = cumsum(n_elements);
if (opt.verbose)
    figure; hist(y, x);
    set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    set(gca, 'XTick', x); set(gca, 'xlim', [-0.5 options.bw(1)+0.5]);
    xlabel('Distance in pixels'); ylabel('Count of occurrences');
    set(gca, 'YTick', 0:50:2*c_elements(end));
    title('Histogram of max errors');
    if (opt.save)
        saveas(gcf, './img/histogram_max.png');
        saveas(gcf, './img/histogram_max.fig');
    end;
    
    figure; bar(x, c_elements / max(c_elements));
    xlabel('RMS [px]'); ylabel('Count of occurrences [%]');
    set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    set(gca, 'XTick', x); set(gca, 'xlim', [-0.5 options.bw(1)+0.5]);
    set(gca, 'YTick', 0:0.05:1);
    title('Cumulative historgam of max errors');
    if (opt.save)
        saveas(gcf, './img/cumhist_max.png');
        saveas(gcf, './img/cumhist_max.fig');
    end;
end;
% save histograms
save('./results/histograms.mat', 'hist_mean', 'hist_max');

%% Show mean and max error for each point from ground truth

mean_nose = mean(err_nose);
mean_canthus_rl = mean(err_canthus_rl);
mean_canthus_lr = mean(err_canthus_lr);
mean_mouth_corner_r = mean(err_mouth_corner_r);
mean_mouth_corner_l = mean(err_mouth_corner_l);
mean_canthus_rr = mean(err_canthus_rr);
mean_canthus_ll = mean(err_canthus_ll);

max_nose = max(err_nose);
max_canthus_rl = max(err_canthus_rl);
max_canthus_lr = max(err_canthus_lr);
max_mouth_corner_r = max(err_mouth_corner_r);
max_mouth_corner_l = max(err_mouth_corner_l);
max_canthus_rr = max(err_canthus_rr);
max_canthus_ll = max(err_canthus_ll);

fprintf('\n\n__________________________________________________________\n');
fprintf('\nResults:\n');
fprintf('Mean error on nose: \t\t\t\t%.4f\t\n', mean_nose);
fprintf('Mean error on mean_canthus_rl: \t\t%.4f\t\n', mean_canthus_rl);
fprintf('Mean error on mean_canthus_lr: \t\t%.4f\t\n', mean_canthus_lr);
fprintf('Mean error on mean_mouth_corner_r: \t%.4f\t\n', mean_mouth_corner_r);
fprintf('Mean error on mean_mouth_corner_l: \t%.4f\t\n', mean_mouth_corner_l);
fprintf('Mean error on mean_canthus_rr: \t\t%.4f\t\n', mean_canthus_rr);
fprintf('Mean error on mean_canthus_ll: \t\t%.4f\t\n', mean_canthus_ll);
fprintf('__________________________________________________________\n');

fprintf('\nR_mean_tst: \t\t\t\t\t%.5f\n', Rtst);
fprintf('R_max_tst: \t\t\t\t\t\t%.5f\n', sum(err_maxdist)/length(err_maxdist));

fprintf('__________________________________________________________\n');

%%

xSC = 0:ceil(max(L));
ySC = round(L);
hist_meanSC = hist(ySC, xSC);
% cumulative histogram of mean errors
n_elementsSC = histc(ySC, xSC); 
c_elementsSC = cumsum(n_elementsSC);
% max
xSCmax = 0:ceil(max(L));
ySCmax = round(err_maxdist);
hist_maxSC = hist(ySCmax, xSCmax);
% cumulative histogram of max errors
n_elementsSCmax = histc(ySCmax, xSCmax); 
c_elementsSCmax = cumsum(n_elementsSCmax);

fig = figure; 
y1 = c_elementsSC / max(c_elementsSC) * 100;
plot(xSC(1:end), y1(1:end), 'b', 'LineWidth', 3); hold on; grid on;
set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
set(gca, 'XTick', xSC); set(gca, 'YTick', 0:5:100, 'FontSize', 15); set(gca, 'ylim', [0 105]);
l1 = legend('corners dataset, normalized trn loss', 'Location', 'SouthEast');
xlabel(['relative error [%]'], 'FontSize', 20); 
ylabel('Count of occurrences [%]', 'FontSize', 20);
title('Cumulative historgam of mean errors', 'FontSize', 20);

% r = 781; c = 731;
r = 801; c = 801;
% Sets position and size of figure on the screen
set(fig, 'Units', 'pixels', 'position', [0 0 c r] ); 
% Sets axes to fill the figure space
% set(gca, 'Units', 'pixels', 'position', [0 0 c+1 r+1 ]);
% Sets print properties; Looks like 1 pixel = (3/4)th of a point
set(fig, 'paperunits', 'points', 'papersize', [fix((c-1)*(3/4))+1 fix((r-1)*(3/4))+1]);
set(fig, 'paperunits', 'normalized', 'paperposition', [0 0 1 1]);
print( fig, sprintf('-r%d', ceil(72*(4/3))), '-dpng', './img/cumhist-mean.png');

fig = figure; 
y2 = c_elementsSCmax / max(c_elementsSCmax) * 100; 
plot(xSCmax(1:end), y2(1:end), 'b', 'LineWidth', 3); hold on; grid on;
set(gcf, 'OuterPosition', [scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
set(gca, 'XTick', xSCmax); set(gca, 'YTick', 0:5:100, 'FontSize', 15);
set(gca, 'ylim', [0 105]);
l1 = legend('corners dataset, normalized trn loss', 'Location', 'SouthEast');
xlabel(['relative error [%]'], 'FontSize', 20); 
ylabel('Count of occurrences [%]', 'FontSize', 20);
title('Cumulative historgam of max errors', 'FontSize', 20);

r = 801; c = 801;
% Sets position and size of figure on the screen
set(fig, 'Units', 'pixels', 'position', [0 0 c r] ); 
% Sets axes to fill the figure space
% set(gca, 'Units', 'pixels', 'position', [0 0 c+1 r+1 ]);
% Sets print properties; Looks like 1 pixel = (3/4)th of a point
set(fig, 'paperunits', 'points', 'papersize', [fix((c-1)*(3/4))+1 fix((r-1)*(3/4))+1]);
set(fig, 'paperunits', 'normalized', 'paperposition', [0 0 1 1]);
print( fig, sprintf('-r%d', ceil(72*(4/3))), '-dpng', './img/cumhist-max.png');

%% Timestamp

fprintf(1,'Finished on %s\n\n', datestr(now));