%% TMP

addpath('./libqp/matlab/');

%%
% Example: Training two-class SVM classifier with L2-regularized bias term.

% load training and testing data
load('riply_dataset','trn_X','trn_y','tst_X','tst_y');

% ensure that the labels are +1/-1
trn_y(find(trn_y~=1)) = -1;
tst_y(find(tst_y~=1)) = -1;

% input arguments
lambda_bmrm = 1e-2;        % regularization parameter 
options.verb = 1;     % will display progress info

nExamples = length(trn_y);
NUM_THREADS = 2;
% TMP
try
   matlabpool('close');
catch ME2
   idSegLast = regexp(ME2.identifier, '(?<=:)\w+$', 'match');
   if (strcmp(idSegLast, 'OpenConnection'))
       fprintf('Cannot close matlabpool. Skipping...\n');
   end;
end;
fprintf('Starting MATLAB Pool...\n');
try
   matlabpool('open', 'local', NUM_THREADS);
catch ME1
   idSegLast = regexp(ME1.identifier, '(?<=:)\w+$', 'match');
   if (strcmp(idSegLast, 'OpenConnection'))
       fprintf('Already connected to matlab pool. Skipping...\n');
   end;
end
%
idx = randperm(nExamples);
data = cell(NUM_THREADS,1);
to = 0;
for i=1:NUM_THREADS    
    from=to+1;
    to = round(i*nExamples/NUM_THREADS);
    data{i}.X0 = 1;
    data{i}.X = trn_X(:,idx(from:to));
    data{i}.y = trn_y(idx(from:to));
end

% lambda = lambda_bmrm*NUM_THREADS;
lambda = lambda_bmrm;

% call BMRM solver
%[W,stat]= parbmrm(data,@risk_svm,lambda,options);
[W,stat]= parbmrm(data,@risk_svm_normalized,lambda,options);

if data{1}.X0 == 0
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
figure;
ppatterns(trn_X,trn_y);
pline(W,W0);
