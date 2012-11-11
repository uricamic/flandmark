function [ aabb ] = makeAABB( bbox )
%MAKEAABB creates AABB - matrix [2 x 4] [x1 ... x4; y1 ... y4] each column
%is one corner point on AABB
%   Detailed explanation goes here
% 
% INPUT:
%   bbox    ...     bbox vector [1 x 4] form [xmin, ymin, xmax, ymax]
% 
% OUTPUT:
%   aabb    ...     AABB matrix [2 x 4], form [x1 .. x4; y1 .. y4]
% 
% 09-07-10 Michal Uricar

     aabb = [bbox(1) bbox(3) bbox(3) bbox(1);
             bbox(2) bbox(2) bbox(4) bbox(4) ];

end

