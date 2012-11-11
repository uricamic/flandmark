function [ bad_idx ] = prepLBPconf( config, options, type, annotation_struct, images_fname, gt_fname )
%PREPLBPCONF Prepares configuration for LBP features function call. 
%   Detailed explanation goes here
% 
% INPUT:
%   config       ...     structure config (see experiment02.m)
%   options      ...     structure options (bw, bw_margin, components, M, order, ...)
%   type         ...     string with one of ['TRN', 'TST', 'VAL'] value
%   datapath     ...     path to dataset (TRN, VAL, TST)
%   images_fname ...     full path filename where matrix with images will be stored
%   gt_fname     ...     full path filename where GT will be stored
% 
% 09-08-10 Michal Uricar
% 12-07-11 Michal Uricar, corners dataset modification

    bad_idx = [];

    if (~strcmp(type, 'TRN') && ~strcmp(type, 'VAL') && ~strcmp(type, 'TST'))
        error('??? Invalid type of dataset. Must be one of [''TRN'', ''VAL'', ''TST'']');
    end;
    
    fprintf(['Prepare ' type ' images...\n']);

%     load(datapath);
    switch (config.nImages)
        case 'all'
            nImages = annotation_struct.N;
        case 'half'
            nImages = round(annotation_struct.N/2);
        otherwise
            nImages = str2num(config.nImages);
    end;
    
    gt = cell(nImages, 1);
    
    fprintf([type ': Creating array of images...\n']);
    imSize = [options.bw(2), options.bw(1)];
    Images = uint8(zeros(options.bw(1)*options.bw(2), nImages));
    for i = 1 : nImages
        [I, Annotation] = getImageFrame(options, i, annotation_struct);
        
        if (isempty(Annotation))
            bad_idx = [bad_idx i];
            continue;
        end;
        
        Points = prepareS0gt(Annotation.P, options);    % transform nose to the center of face
        Points = Points(:, options.comselect);          % extract relevant points only (name list in options.compnames)
        Points(:, options.M) = Annotation.P(:, 10);     % copy back original nose position

        fprintf([type ': %.2f%% Processing image %s.\n'], i*100/nImages, Annotation.image.filename);
        
        Images(:, i) = I(:);
        gt{i} = Points;
        %gt{i} = TestPoints;
    end;
    
    %% Save array of Images
    fprintf([ type ': Saving array of images to file %s...\n'], images_fname);
    save(images_fname, 'Images', 'imSize');
    fprintf([ type ': Saving ground truth points to file %s...\n'], gt_fname);
    save(gt_fname, 'gt');
end
