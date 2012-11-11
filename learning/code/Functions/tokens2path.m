function [ spath ] = tokens2path( tokens )
%TOKENS2PATH Summary of this function goes here
%   Detailed explanation goes here

    spath = tokens{1};

    for i = 2 : numel(tokens)
        spath = [spath '/' tokens{i}];
    end;

end
