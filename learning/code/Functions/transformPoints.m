function [ Pt ] = transformPoints( P, bb, bw )
%TRANSFORMPOINTS Transforms coordinates of detected landmarks into the
%normalized image frame.
%   Detailed explanation goes here
% 
% INPUT:
%   P       ...     [2 x M (double)] Detected landmarks coordinates in 
%                   normalized image frame. M is the number of components
%   bb      ...     [1 x 4 (double)] bbox of normalized image frame in
%                   original image (i.e. the extended bbox of face
%                   detection)
%   bw      ...     [2 x 1 (double)] base window size
% 
% OUTPUT:
%   Pt      ...     [2 x M (double)] Tranformed coordinates of detected
%                   landmarks 
% 
% 05-08-11 Michal Uricar

    Origin = [bb(1); bb(2)];
    
    % compute the scale factor
    scalefac = [ (bb(3)-bb(1)+1)/bw(1);
                 (bb(4)-bb(2)+1)/bw(2) ];
             
    M = size(P, 2);
    
    Pt = P .* repmat(scalefac, 1, M) + repmat(Origin, 1, M);

end

