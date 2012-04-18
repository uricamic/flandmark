% This script demonstrates SVMOCAS_LBP solver for training two-class 
% classifiers of gray scale images. The images are described by a pyramid
% of Local Binary Patterns (LBP) features. In order to make the classifier
% translation invariant, SVMOCAS_LBP allows to generate virtual examples
% by shifting a base window along x- and y-axis of the input images.
% SVMOCAS_LBP implements the COFFIN framework, i.e. the virtual examples 
% are computed on demand which leads to substantial memory savings. 
%
% For more details see:
%  S. Sonnenburg, V. Franc.  COFFIN: A Computational Framework for Linear SVMs.
%  In Proc. of ICML 2010. Haifa. 
%
% Image database provided by courtesy of 
%  Eyedea Recognition LtD. (http://www.eyedea.cz)
%

load('./data/gender_images.mat',...
     'trn_male_images','trn_female_images',...
     'tst_male_images','tst_female_images','IMAGE_SIZE');

%% Parameter of LBP feature representation
LBP_PYRAMID = 4;       % height of LBP pyramid; 1 means that LBP filters are
                       % applied only on the original image otherwise
                       % the image is LBP_PYRAMID-times downscaled
 
%% Definition of base window on which the classifier is applied 
WINDOW_SIZE = [60 40];     % [height width]
WINDOW_TOP_LEFT_COL = 20;
WINDOW_TOP_LEFT_ROW = 15;

%% Parameters SVM solver
X0 = 1;                 % added constant feature
C = 0.001;              % SVM C
SVMOCAS_BUFSIZE = 100;  % number of cutting planes to store; single CP requires ~ 0.7MB
SVMOCAS_TOLREL = 0.01;  % precision of SVM solution
SVMOCAS_METHOD = 1;     % 1.. OCAS, 0..BMRM

%% Virtual examples
% The virtual examples are created by shifting the base window in x and
% y-axis and mirroring the image along its vertical axis.

switch input('use virtual examples (0-no, 1-yes)?:')
    case 0
        % testing error 20%
        USE_VIRTUAL_EXAMPLES = 0; 

    case 1
        % testing error 14%
        USE_VIRTUAL_EXAMPLES = 1; 
        SHIFT_X = [-1 0 1];
        SHIFT_Y = [-1 0 1];
        MIRROR = [0];
        
end

%% prepare training examples 
nTrnMales = size(trn_male_images,2);
nTrnFemales = size(trn_male_images,2);

if ~USE_VIRTUAL_EXAMPLES
    BOX = [WINDOW_TOP_LEFT_COL;WINDOW_TOP_LEFT_ROW;0];
    wins = [ [1:nTrnMales nTrnMales+[1:nTrnFemales]]; repmat(BOX,1,nTrnFemales+nTrnMales)];
    labels = [ones(1,nTrnMales) -ones(1,nTrnFemales)];    
        
else               
    % define virtual examples by shifting the base window
    BOX = zeros(3,length(MIRROR)*length(SHIFT_Y)*length(SHIFT_X),'uint32');    
    cnt = 0;
    for u=SHIFT_X
        for v=SHIFT_Y
           for mirror=MIRROR
               cnt = cnt + 1;
               BOX(:,cnt) = uint32([WINDOW_TOP_LEFT_COL+u; WINDOW_TOP_LEFT_ROW+v;mirror]);
           end
        end
    end       

    % prepare virtual examples for each training image
    wins = zeros(4,size(BOX,2)*(nTrnMales+nTrnFemales),'uint32');
    labels = zeros(size(BOX,2)*(nTrnMales+nTrnFemales),1);
    cnt = 0;
    for i=1:(nTrnMales+nTrnFemales)
        wins(1,cnt+1:cnt+size(BOX,2)) = uint32(i);
        wins(2:4,cnt+1:cnt+size(BOX,2)) = BOX;
        if i <= nTrnMales
            labels(cnt+1:cnt+size(BOX,2)) = 1;
        else
            labels(cnt+1:cnt+size(BOX,2)) = -1;
        end
        cnt = cnt + size(BOX,2);
    end
    
end

%% train SVM linear classifier 
[W,W0,stat] = svmocas_lbp([trn_male_images trn_female_images], IMAGE_SIZE,...
                          uint32(wins), WINDOW_SIZE, LBP_PYRAMID, X0, labels, C, ...
                          1, SVMOCAS_TOLREL,0,0,SVMOCAS_BUFSIZE);
                          
%% compute features for testing images and evaluate the tained classifier
BOX = [WINDOW_TOP_LEFT_COL;WINDOW_TOP_LEFT_ROW;0];

nTstMale = size(tst_male_images,2);
tst_labels = ones(1,nTstMale);
tst_feat = lbppyr_features(tst_male_images,IMAGE_SIZE, ...
                          uint32([1:nTstMale;repmat(BOX,1,nTstMale)]), ...
                          WINDOW_SIZE, LBP_PYRAMID,0);

nTstFemale = size(tst_female_images,2);
tst_labels = [tst_labels -ones(1,nTstFemale)];
tst_feat = [tst_feat lbppyr_features(tst_female_images,IMAGE_SIZE, ...
                          uint32([1:nTstFemale;repmat(BOX,1,nTstFemale)]), ...
                          WINDOW_SIZE, LBP_PYRAMID,0)];
                                                
ypred = sign(W'*double(tst_feat) + W0);
tst_error = sum(ypred(:) ~= tst_labels(:))/length(tst_labels);
fprintf('Testing error: %f%%\n', tst_error*100);
                      