function [ face_image, bb ] = deprecated_getNormalizedFrame( I, bbox, options )
%GETNORMALIZEDFRAME crops image I to bbox extended by options.bw_margin and
%scales it to options.bw
%   Detailed explanation goes here
% 
% INPUT:    
%   I       ...     input image
%   bbox    ...     bounding box provided by face detector
%   options ....    structure with parameters bw and bw_margin
% 
% OUTPUT
%   face_image  ...     normalized image frame (input for function detector)
% 
% 20-10-10 Michal Uricar, fixed d computation and imcrop rect
% 18-08-10 Michal Uricar

    % extend bbox by bw_margin
    d = [bbox(3)-bbox(1)+1; bbox(4)-bbox(2)+1];                         % get bbox dimensions
    C = [(bbox(3)+bbox(1))/2; (bbox(4)+bbox(2))/2];                     % get center of bbox
    nd = d .* options.bw_margin./100 + d;                               % get extended dimensions (by bw_margin)
    bb = [ C(1)-nd(1)/2,  C(2)-nd(2)/2,  C(1)+nd(1)/2,  C(2)+nd(2)/2 ]; % create extended bbox

    % crop and resize image to normalized frame
    face_image = imcrop(I, [bb(1), bb(2), (bb(3)-bb(1)+1), (bb(4)-bb(2)+1)]);
    face_image = imresize(face_image, [options.bw(2) options.bw(1)]); 
end
