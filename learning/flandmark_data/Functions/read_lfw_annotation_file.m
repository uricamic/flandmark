function [ annotaion_struct ] = read_lfw_annotation_file( filename )
%READ_LFW_ANNOTATION_FILE Summary of this function goes here
%   Detailed explanation goes here

    fid = fopen(filename, 'r');
    C = textscan(fid, '%s %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
    fclose(fid);

    annotaion_struct = [];
    annotaion_struct.names = C{1};
    annotaion_struct.bbox = [C{2} C{3} C{4} C{5}];
    annotaion_struct.eye_r = [C{6} C{7}];
    annotaion_struct.eye_l = [C{8} C{9}];
    annotaion_struct.canthus_rr = [C{10} C{11}];
    annotaion_struct.canthus_rl = [C{12} C{13}];
    annotaion_struct.canthus_lr = [C{14} C{15}];
    annotaion_struct.canthus_ll = [C{16} C{17}];
    annotaion_struct.mouth = [C{18} C{19}];
    annotaion_struct.mouth_corner_r = [C{20} C{21}];
    annotaion_struct.mouth_corner_l = [C{22} C{23}];
    annotaion_struct.nose = [C{24} C{25}];
    annotaion_struct.N  = numel(annotaion_struct.names);
    
end

