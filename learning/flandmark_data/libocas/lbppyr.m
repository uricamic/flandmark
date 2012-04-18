function F = lbppyr(I, P)
% LBPPYR computes LBP features on scale-pyramid of input image.
%
% Synopsis:
%  F = lbppyr(I, P)
%
% Description:
%  Local Binary Pattern (LBP) is 1 Byte number which encodes grey-scale
%  intensities of a given 3x3 window.
% 
%  For P = 1, LBPPYR computes LBP features for all 3x3 windows
%  in the input image I. The LBP features are stucked to a column 
%  vector F.
%
%  For P > 1, the input image I is down-scaled (P-1)-times to create 
%  a scale-pyramid of height P. Each subsequent level of the pyramid is
%  obtained by downscaling the original image by 2. Having the pyramid,
%  the LBP features are computed on each leavel of the pyramid and 
%  they are stucked to a column vector F.
%
%  For more details on LBP see:
%   Ojala, et. al: Multiresolution gray-scale and rotation invariant 
%   texture classification with local binary patterns. IEEE PAMI, 
%   24(7):971-987,2002.
%
% Input:
%  I [H x W (uint8)] Input image.
%  lbp_pyramid [1 x 1] Height of the LBP pyramids.
%  
% Output:
%  F [nDim x 1] LBP features stucked to a column vector.
%
% Example:
%  I = imread('./data/lena.jpg');
%  F = lbppyr(I,4);
%
% Matlab code equivalent to LBPPYR_MEX.C:
% 
%  [h,w] = size(I);
%  I = uint32(I);
%  F = [];
%  k = 1; 
%  while k <= P & min(w,h) >= 3
%     k = k + 1;
%     tmp = lbpfilter(double(I));
%     tmp = tmp(2:end-1,2:end-1);
%     F = [F; tmp(:)];    
%     if mod(w,2) == 1
%        I = I(:,1:w-1);
%        w= w - 1;
%     end
%     if mod(h,2) == 1
%        I = I(1:h-1,:);
%        h = h - 1;
%     end
%     I = I(1:2:h,1:2:w)+I(1:2:h,2:2:w)+I(2:2:h,1:2:w)+I(2:2:h,2:2:w);
%     [h,w] = size(I);
%  end
%