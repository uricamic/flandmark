function [ L ] = computeL( options, GT, kappa)
%COMPUTEL computes values of penalty function L(S, S*) for each component
%   Detailed explanation goes here
% 
% INPUT:
%   options     ...     options (S, components, M)
%   GT          ...     matrix [ 2 x M ] ground truth position of landmarks
% 
% OUTPUT:
%   L           ...     cell array [ M x 1 ], values of penalty function L(S, S*)  for each component
% 
% 30-07-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset (checked-only)
% 29-07-11 Michal Uricar, loss function with kappa normalization
    
    L = cell(options.M, 1);

    for comp = 1 : options.M
        % size of component's bbox
        siz = [options.S(4, comp) - options.S(2, comp) + 1; options.S(3, comp) - options.S(1, comp) + 1];
        % offset (top left corner of component's bbox)
        offset = options.S(1:2, comp);
        
        cor_sqr = zeros(2, siz(2)*siz(1));
        
        for i = 1 : siz(2)
            % x coordinates of (S_i - S_i*)^2
            cor_sqr(1, (i-1)*siz(1)+1:(i*siz(1))) = (i-1+offset(1) - GT(1, comp)).^2;
            % y coordinates of (S_i - S_i*)^2
            cor_sqr(2, (i-1)*siz(1)+1:(i*siz(1))) = (offset(2) - GT(2, comp):offset(2)+siz(1)-1 - GT(2, comp)).^2;
        end;
        % L(S, S*) = 1/M \sum_{i=0}^{M-1} || S_i - S_i* ||
        if (nargin < 3)
            L{comp} = 1/options.M .* sqrt( cor_sqr(1,:) + cor_sqr(2,:) );
        else 
            L{comp} = 100 * kappa * 1/options.M .* sqrt( cor_sqr(1,:) + cor_sqr(2,:) );
        end;
    end;

end
