function [ gt ] = prepareS0gt( points, options )
%PREPARES0GT Summary of this function goes here
%   
%   make line between eyes, find center, make line perpendicular to line
%   connecting eyes, do perpendicular projection of mouth to this line and
%   then find point in 1/2 of length of segment connecting projected mouth
%   end center of eyes. This point is new GT position of S0 component.
% 
% 11-07-11 Michal Uricar, corners dataset modification

    gt = points;

    eyer = points(:, 1);
    eyel = points(:, 2);
    m1 = points(:, 7);
    
    % centeroid of eyes
    c1 = (eyer + eyel)/2;
    % line between eyel and eyer
    e = crossprod([eyer; 1], [eyel; 1]);
    % shift centeroid of eyes (c1) along line perpendicular to e
    c2 = c1 + e(1:2);
    % create line perpendicular to e trough c1 (and c2 also)
    p = crossprod([c1; 1], [c2; 1]);
    % shift mouth along normal vector of p
    m2 = m1 + p(1:2);
    % create line m passing trough m1 and m2
    m = crossprod([m1; 1], [m2; 1]);
    % compute intersection of lines m and p and normalize
    m3 = crossprod(m, p); m3 = m3./m3(3);
    
    % S0 is centeroid of c1 and m3
    gt(:, 10) = round( (c1 + m3(1:2))/2 );
    
    % compute distance of mouth to line p
    %d = (p(1:2)'*mouth + p(3))/norm(p(1:2));
end

