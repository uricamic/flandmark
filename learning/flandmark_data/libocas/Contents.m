% OCAS solver for training linear SVM classifiers from large-scale data
%
% Copyright (C) 2008,2009,2010 
% Vojtech Franc, xfrancv@cmp.felk.cvut.cz
% Soeren Sonnenburg, soeren.sonnenburg@tu-berlin.de
%
% SVM solvers:
%   svmocas             Train linear binary (two-class) SVM classifier using 
%                       OCAS solver.
%   svmocas_lbp         Train linear SVM classifier for images represented by 
%                       LBP features. 
%   msvmocas            Train multi-class linear SVM classifier using OCAS solver.
%
% Auxciliary functions:
%   svmlight_linclass   Classify examples in SVM^light file by linear rule.
%   svmocas_parseout    Parsing of text output of SVMOCAS solver.
%
% Examples
%   libocas_test        This script tests functionality of SVMOCAS and 
%                       MSVMOCAS solvers.
%   svmocas_lbp_example Example on using SVMOCAS_LBP from training translation
%                       invariant image classifiers.
%
% LBP features for image classification:
%   lbpfilter           Computes LBP for each 3x3 subwindow in input image.
%   lbppyr              Computes LBP features on scale-pyramid of input image.
%   lbppyr_features     Computes pyramid of LBP features for each defined 
%                       window in input images. 
%