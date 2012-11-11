function [ Psi ] = buildPsi2( options, psi, configuration )
%BUILDPSI functions creates vector \Psi(I, S) for given configuration
%   Detailed explanation goes here
% 
% INPUT
%   options         ...     structure holds info about components (size, count) and S (statistics of component position)
%   psi             ...     cell Mx1 contains matrices \Psi^q(I, S_i)
%   configuration   ...     matrix 2xM [x0, ..., x_{M-1}; y_0, ..., y_{M-1}] with positions of components
% 
% OUTPUT
%   Psi     ...     vector [\Psi_0^q(I, S_0); ... ; \Psi_{M-1}^q(I, S_{M-1}); g(S_0xS_1); ... ; g(S_0xS_{M-1})];
% 
% 10-08-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset modificatoin

    indices = zeros(1, options.M);
    
    Psi = [];

    for i = 1 : options.M
        psi_si = psi{i};
        P = configuration(:, i) - [options.S(1, i) - 1; options.S(2, i) - 1];
        siz = [options.S(4, i) - options.S(2, i) + 1; options.S(3, i) - options.S(1, i) + 1];
        indices(i) = (P(1) - 1)*siz(1) + P(2);
        Psi = [Psi; psi_si(:, indices(i))];
    end;
    
    % g
    %for i = 2 : options.M - 2
    for i = 2 : options.M - 3
        dxdy = configuration(:, i) - configuration(:, 1);
        g = [dxdy; dxdy.^2];
        Psi = [Psi; g];
    end;
    % g5
    dxdy = configuration(:, 6) - configuration(:, 2);
    g = [dxdy; dxdy.^2];
    Psi = [Psi; g];
    % g6
    dxdy = configuration(:, 7) - configuration(:, 3);
    g = [dxdy; dxdy.^2];
    Psi = [Psi; g];
    % g7
    dxdy = configuration(:, 8) - configuration(:, 1);
    g = [dxdy; dxdy.^2];
    Psi = [Psi; g];
    
end
