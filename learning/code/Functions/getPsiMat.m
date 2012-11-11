function [ M, N ] = getPsiMat( data, image )
%GETPSIMAT gets \Psi matrix (LBP features) for all components
%   Detailed explanation goes here
% 
% Synopsis:
%  [ M ] = getPsiMat(data, image)
% 
% INPUT:
%   data    ...     struct, containes filed lbp (cell array - holds
%                   neccessary parameters for lbppyr_features_sparse
%                   function), imSize (size of image), Images (matrix
%                   [imSize(1)*imSize(2) x nImages]), nImages (count of
%                   images), structure options (with field M - count of
%                   components)
%   image   ...     idx of image on which LBP features are computed
% 
% OUTPUT:
%   M       ...     cell array [1 x data.options.M] - \Psi matrices with
%                   LBP features
%   N       ...     cell array [1 x data.options.M] - \Psi matrices with
%                   indices of nonzero elements of LBP features
% 
% 05-08-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset (checked-only)

    M = cell(1, data.options.M);
    N = cell(1, data.options.M);

    for j = 1 : data.options.M;
        [M{j} N{j}] = lbp2sparse(data.Images(:, image), data.imSize, data.lbp{j}.wins, data.lbp{j}.winSize, data.lbp{j}.hop, 0);
    end;
    
end

