function [ detection, Q, G ] = detector( face_image, model )
%DETECTOR Given normalized image frame face_image and model structure this
%function detects positions of significant points on human face and returns
%its positions.
% 
% INPUT: 
%   face_image      ...     normalized image frame (see function getNormalizedFrame)
%   model           ...     structure holding parameters of model (see
%                           exampleUsageDetector.m for details how to create it) 
% 
% OUTPUT:
%   detection       ...     matrix 2xM format [x_0 ... x_{M-1}; 
%                                              y_0 ... y_{M-1}]
%                           detected positions of significant points
% 
% 18-08-10 Michal Uricar
% 08-08-11 Michal Uricar, argmax_mex modification

    model.data.Images = face_image(:);

%     Psi = model.fgetPsiMat(model.data, 1);
%     Psi = model.fgetPsiMat(model.data);

%     [ q, g ] = model.fgetQG(model.data.options, model.W, Psi, model.data.mapTable);
%     detection = model.fargmax(model.data.options, q, g);

    [Psi, PsiSparse] = model.fgetPsiMat(model.data, 1);

    if (nargout > 1)
        [ q, g ] = model.fgetQG(model.data.options, model.W, Psi, model.data.mapTable);
        [detection, Q, G] = model.fargmax(model.data.options, q, g);
    else
        %detection = model.fargmax(model.data.options, q, g);
        detection = model.argmax_mex(model.data.options, model.W, model.data.mapTable, PsiSparse);
    end
    
end

