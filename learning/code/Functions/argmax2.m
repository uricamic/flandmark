function [ S_max ] = argmax2( options, q, g )
%ARGMAX Function S = argmax_{S \in \mathcal{S}} f(I, S)
%   Detailed explanation goes here
% 
% INPUT:
%   options        ...     options (components, S, M)
%   q              ...     cell [1xM]  values of function q_i(I, S_i) = <w_i^q, \Psi_i^q(I, S_i)>
%   g              ...     cell [1xM-1] values of function g_i(S_0, S_i) = <w_i^g, \Psi_i^g(S_0, S_i)>
%   options.PsiGS0 ...     cell matrix [|S0| x M-1] of displacements
%   options.PsiGS1 ...     cell matrix [|S1| x 1] of displacements
%   options.PsiGS2 ...     cell matrix [|S2| x 1] of displacements
% 
% OUTPUT:
%   S_max          ...     matrix 2xM [x_0, ... , x_{M-1}; y_0, ..., y_{M-1}] positions of landmarks
%
% 10-08-10 Michal Uricar
% 17-10-10 Michal Uricar, added computation of Q and G for output
% 11-07-11 Michal Uricar, corners dataset

    S_max = zeros(2, options.M);
    q0 = q{1}; q0_length = length(q0);
    s0 = zeros(options.M, q0_length);
    
    % store maximum and index of s5 for all positions of s1
    q1_length = length(q{2});
    s1 = zeros(2, q1_length);
    for i = 1 : q1_length
        g5 = g{5}'*options.PsiGS1{i, 1};
        [s1(1, i), s1(2, i)] = max( q{6}' + g5 );
    end;
    s1(1, :) = s1(1, :) + q{2}';
    % store maximum and index of s6 for all positions of s2
    q2_length = length(q{3});
    s2 = zeros(2, q2_length);
    for i = 1 : q2_length
        g6 = g{6}'*options.PsiGS2{i, 1};
        [s2(1, i), s2(2, i)] = max( q{7}' + g6 );
    end;
    s2(1, :) = s2(1, :) + q{3}';
    
    for i = 1 : q0_length
        sums0 = 0;
        q10 = s1(1, :) + g{1}'*options.PsiGS0{i, 1};
        [maxim, s0(2, i)] = max( q10 );
        sums0 = sums0 + maxim;
        s0(6, i) = s1(2, s0(2, i));
        
        q20 = s2(1, :) + g{2}'*options.PsiGS0{i, 2};
        [maxim, s0(3, i)] = max( q20 );
        sums0 = sums0 + maxim;
        s0(7, i) = s2(2, s0(3, i));
        
        q30 = q{4}' + g{3}'*options.PsiGS0{i, 3};
        [maxim, s0(4, i)] = max( q30 );
        sums0 = sums0 + maxim;
        
        q40 = q{5}' + g{4}'*options.PsiGS0{i, 4};
        [maxim, s0(5, i)] = max( q40 );
        sums0 = sums0 + maxim;
        
        q70 = q{8}' + g{7}'*options.PsiGS0{i, 5};
        [maxim, s0(8, i)] = max( q70 );
        sums0 = sums0 + maxim;

        s0(1, i) = sums0;
    end;
    s0(1, :) = s0(1, :) + q0';

    [maxim, idx] = max(s0(1, :));
    indices = [idx; s0(2:end, idx)];
    
    % transform indices to 2D coordinates
    for i = 1 : options.M
        siz = [options.S(4, i) - options.S(2, i) + 1; options.S(3, i) - options.S(1, i) + 1];
        ind_i = mod(indices(i) - 1, siz(1)); ind_j = floor((indices(i) - 1) / siz(1));
        S_max(:, i) = [ind_j + options.S(1, i); ind_i + options.S(2, i)];
    end;    
end
