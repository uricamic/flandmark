function [ mapTable ] = createMapTable2( options, Data )
%CREATEMAPTABLE creates map table to easily extract w_q^i and w_g^i from
%vector W.
%   Detailed explanation goes here
% 
% INPUT:
%   options     ...     structure (must contain field M - count of
%                       components) 
%   Data        ...     cell array - Psi matrices with LBP features for
%                       each component
% 
% OUTPUT:
%   mapTable    ...     matrix [options.M x 4] format [start_q_i, end_q_i,
%                       start_g_i-1, end_g_i-1; ... ]
% 
% 10-08-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset (checked-only)

    mapTable = zeros(options.M, 4);

    from = 1; to = size(Data{1}, 1);
    mapTable(1, 1) = from; mapTable(1, 2) = to;
    
    for i = 2 : options.M
        from = to + 1;
        to = from + size(Data{i}, 1) - 1;
        mapTable(i, 1) = from; mapTable(i, 2) = to;
    end;
    
    mapTable(2:end, 3) = to+1:4:to+(options.M-1)*4;
    mapTable(2:end, 4) = to+4:4:to+(options.M-1)*4;
        
end

