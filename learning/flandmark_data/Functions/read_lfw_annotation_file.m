function [ annotaion_struct ] = read_lfw_annotation_file( filename )
%READ_LFW_ANNOTATION_FILE Summary of this function goes here
%   Detailed explanation goes here
% 
% 
%   20-04-12 Michal Uricar, version for LFW annotation with 7 points
%                           (canthi, mouth corners and nose)

    fid = fopen(filename, 'r');
    C = textscan(fid, '%s %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);
    
    annotaion_struct = [];
    annotaion_struct.names = C{1};
    annotaion_struct.bbox = [C{2} C{3} C{4} C{5}];
    annotaion_struct.canthus_rr = [C{6} C{7}];
    annotaion_struct.canthus_rl = [C{8} C{9}];
    annotaion_struct.canthus_lr = [C{10} C{11}];
    annotaion_struct.canthus_ll = [C{12} C{13}];
    annotaion_struct.mouth_corner_r = [C{14} C{15}];
    annotaion_struct.mouth_corner_l = [C{16} C{17}];
    annotaion_struct.nose = [C{18} C{19}];
    annotaion_struct.N  = numel(annotaion_struct.names);
   
    % compute coordinates of center of eyes and center of mouth
    annotaion_struct.eye_r = [(C{6}+C{8})/2 (C{7}+C{9})/2];
    annotaion_struct.eye_l = [(C{10}+C{12})/2 (C{11}+C{13})/2];
    annotaion_struct.mouth = [(C{14}+C{16})/2 (C{15}+C{17})/2];
    
end
