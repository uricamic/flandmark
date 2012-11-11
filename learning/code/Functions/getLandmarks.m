function [ Points, Names ] = getLandmarks( face, index )
%GETLANDMARKS returns position of landmarks from structure face in matrix
%form, also returns array of strings - names of landmarks
%   Detailed explanation goes here
% 
% INPUT
%   face    ...     structure with face XML annotation (containes positions of landmarks)
%   index   ...     index where landmarks in structure face actualy begins
% 
% OUTPUT
%   Points  ...     Matrix 2xM of landmarks positions [x_0, ..., x_{M-1}; y_0, ..., y_{M-1}]
%   Names   ...     array of cells, each cell holds string - name of component
% 
% 14-07-10 Michal Uricar
% 21-03-12 Michal Uricar, deprecated
        
    tmp = struct2cell(face);
    
    Points = zeros(2, numel(tmp) - index);
    
    for i = index : numel(tmp)
        Points(:, i - index + 1) = cell2mat(tmp(i, :));
    end;

    Names = fieldnames(face);
    Names(1:index-1) = [];
    
end

