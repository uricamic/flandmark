addpath('./libqp/matlab/');

% Example: Training two-class SVM classifier with L2-regularized bias term.

% load training and testing data
load('riply_dataset','trn_X','trn_y','tst_X','tst_y');

% ensure that the labels are +1/-1
trn_y(find(trn_y~=1)) = -1;
tst_y(find(tst_y~=1)) = -1;

% input arguments
% lambda = 1e-2;        % regularization parameter 
lambda = 10^-0;        % regularization parameter 
options.verb = 1;     % will display progress info

% prepare training examples for BMRM
data.X0 = 1;        % added constant feature -> biased classification rule
data.X = trn_X;     % inputs
data.y = trn_y;     % labels

% call BMRM solver
[W,stat]= bmrm(data,@risk_svm,lambda,options);

if data.X0 == 0
    W0 = 0;
else   
    W0 = W(end);
    W = W(1:end-1);
end

% compute training and testing error
ypred = sign(W'*trn_X + W0);
trn_error = sum(ypred(:)~=trn_y(:))/length(trn_y)

ypred = sign(W'*tst_X + W0);
tst_error = sum(ypred(:)~=tst_y(:))/length(tst_y)

% visualize bounary of the trained liner classification rule
%figure;
%ppatterns(trn_X,trn_y);
%pline(W,W0);
