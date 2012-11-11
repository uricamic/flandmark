function [ q, g ] = getQG( options, W, Data, mapTable, L )
%GETQG function cuts long vector W into corresponding parts w_i^q and
%w_i^g
%   Detailed explanation goes here
% 
% INPUT:
%   options     ...     options (components, S, M)
%   W           ...     long vector W [\Psi_0^q(I, S_0); ... ; \Psi_{M-1}^q(I, S_{M-1}); g(S_0xS_1); ... ; g(S_0xS_{M-1})];
%   Data        ...     matrix of cells (LBP features for each component)
%   mapTable    ...     matrix [ M x 4 ] format [ q0_from, q0_to, 0, 0; q1_from, q1_to, g1_from, g1_to; ... ]
%   L           ...     cell [ M x 1 ] values of penalty function for individual components
% 
% OUTPUT:
%   q           ...     cell [M x 1]  values of function q_i(I, S_i) = <w_i^q, \Psi_i^q(I, S_i)>
%   g           ...     cell [M-1 x 1] values of function g_i(S_0, S_i) = <w_i^g, \Psi_i^g(S_0, S_i)>
%
% 29-07-10 Michal Uricar

    if (nargin < 5)
        L = [];
    end

    q = cell(options.M, 1);
    g = cell(options.M - 1, 1);

    for i = 1 : options.M
        % get q
        if (isempty(L))
            q{i} = (W(mapTable(i, 1):mapTable(i, 2))'*Data{i})';
        else
            q{i} = (W(mapTable(i, 1):mapTable(i, 2))'*Data{i} + L{i})';
        end;
        % get g
        if (i > 1)
            g{i-1} = W(mapTable(i, 3):mapTable(i, 4));
        end;
    end;
    
end
