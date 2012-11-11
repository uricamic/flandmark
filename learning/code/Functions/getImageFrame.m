function [ Iframe, Annotation, I, bbox_orig, bbox, OrigPoints ] = getImageFrame( options, idx, annotation_struct )
%GETIMAGEFRAME return normalized frame and its annotation
%   Given image from db and detected frame of face on it this function
%   expands detected frame by options.bw_margin and than crops image to
%   this frame and rescales it to size of options.bw (base window)
% 
% INPUT
% 
%   options     ...     structure containing parameters - bw (base_window),
%                       bw_margin (base window margin), components (marix with bboxes of components), 
%                       image_path (path to directory containing jpg files of faces)
%   idx         ...     index of image in db
%   image       ...     cell image with filename and  detected frame
% 
% OUTPUT
% 
%   Iframe      ...     detected frame normalized to bw size after
%                       expansion of bw_margin in both dimensions
%   Annotation  ...     structure with annotation in similar format as
%                       original XML annotation
%   I           ...     original Image (rgb)
%   bbox_orig   ...     original bbox of detected face
%   bbox        ...     extended bbox by options.bw_margin
%   OrigPoints  ...     original landmarks (in unnormalized image)
% 
% 02-08-10 Michal Uricar
% 20-10-10 Michal Uricar, normalized image frame get by getNormalizedFrame.m
% 12-07-11 Michal Uricar, corners dataset
% 21-03-12 Michal Uricar, LFW annotation in one file

    M = 10;
    fname = strtok(annotation_struct.names{idx}, '.');
    bbox_orig = double(annotation_struct.bbox(idx, :));
    filename = [options.image_path annotation_struct.names{idx}];
    I = rgb2gray(imread(filename));

    % get normalized image frame
    try
        [Iframe, bbox] = getNormalizedFrame(I, bbox_orig, options);
    catch
        Iframe = [];
        Annotation = [];
        bbox = []; 
        OrigPoints = [];
        return;
    end;

    % compute coordinates of landmarks in normalized frame
    OrigPoints = [annotation_struct.eye_r(idx, :)' annotation_struct.eye_l(idx, :)' ,...
        annotation_struct.canthus_rr(idx, :)' annotation_struct.canthus_rl(idx, :)' ,...
        annotation_struct.canthus_lr(idx, :)' annotation_struct.canthus_ll(idx, :)' ,...
        annotation_struct.mouth(idx, :)' annotation_struct.mouth_corner_r(idx, :)'  ,...
        annotation_struct.mouth_corner_l(idx, :)' annotation_struct.nose(idx, :)'];
    O = [bbox(1); bbox(2)]; 
    scalefac = [ options.bw(1)/(bbox(3)-bbox(1)+1); options.bw(2)/(bbox(4)-bbox(2)+1) ];

    if (size(OrigPoints, 2) ~= M)
        error('??? Count of landmarks does not equal M. Corrupted xml annotation?');
        Annotation = [];
        fprintf('Bad annotation found.\n');
        return;
    end;

    % save annotation
    Points = zeros(2, M);
    for j = 1 : M
        Points(:, j) = (OrigPoints(:, j) - O) .* scalefac;
    end;
    
    if (numel(Points(Points < 0)) > 0)
        Annotation = [];
        fprintf('Bad annotation found.\n');
        return;
    end;
    
    [h w d] = size(Iframe);
    
    Annotation.P = round(Points);
    Annotation.image.filename = [fname '.jpg'];
    Annotation.image.width = w;
    Annotation.image.height = h;
end
