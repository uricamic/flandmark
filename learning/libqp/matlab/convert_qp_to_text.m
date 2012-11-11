% This script converts QP tasks stored in MAT files
% into text file with a simple format.
%
% The text files serve as input to the standalone application examples/qpsplx 
%

Files = {'../data/splx_qp_random100x10','../data/splx_qp_random300x30'};

for i=1:length(Files)
    load([Files{i} '.mat']);
    fid = fopen([Files{i} '.txt'],'w+');

    fprintf(fid,'%d\n', 5);
    fprintf(fid,'H %d %d\n', size(H,1), size(H,2)); 
    fprintf(fid,'%.10f\n', H(:));

    fprintf(fid,'f %d %d\n', size(f,1), size(f,2)); 
    fprintf(fid,'%.10f\n', f(:));

    fprintf(fid,'b %d %d\n', size(b,1), size(b,2)); 
    fprintf(fid,'%.10f\n', b(:));

    fprintf(fid,'I %d %d\n', size(I,1), size(I,2)); 
    fprintf(fid,'%.10f\n', I(:));

    fprintf(fid,'S %d %d\n', size(S,1), size(S,2)); 
    fprintf(fid,'%.10f\n', S(:));

    fclose(fid);
end




