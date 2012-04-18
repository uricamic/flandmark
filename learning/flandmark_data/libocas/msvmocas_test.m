trn_file = 'example4_train.light';
tst_file = 'example4_test.light';

C = 1;
TolRel = 0.01;
TolAbs = 0.00;
QPBound = 0;
BufSize = 2000;
nData = inf;
MaxTime = inf;
verb = 0;

Method = 0; % Cutting Plane Algorithm
fprintf('Training SVM by Cutting Plane Algorithm...');
[cp_W,cp_stat] = msvmocas(trn_file,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb);
fprintf('done\n');

fprintf('Evaluating classifier on testing data...');
[score,true_labels] = svmlight_linclass(tst_file,cp_W,[]);
[dummy,pred_labels] = max(score);
cp_tst_err = sum(pred_labels(:) ~= true_labels(:))/length(true_labels);
fprintf('done\n');


fprintf('Training time: %f[s], #trn_errors: %d, nIter: %d\n',...
        cp_stat.total_time, cp_stat.nTrnErrors, cp_stat.nIter);
fprintf('Objval  primal: %f, dual: %f, gap: %f\n', ...
        cp_stat.Q_P, cp_stat.Q_D,cp_stat.Q_P-cp_stat.Q_D);
fprintf('Testing error: %f %%\n',cp_tst_err*100);
   

Method = 1; % OCAS 
fprintf('\nTraining SVM by OCAS...');
[ocas_W,ocas_stat] = msvmocas(trn_file,C,Method,TolRel,TolAbs,QPBound,BufSize,nData,MaxTime,verb);
fprintf('done\n');

fprintf('Evaluating classifier on testing data...');
[score,true_labels] = svmlight_linclass(tst_file,ocas_W,[]);
[dummy,pred_labels] = max(score);
ocas_tst_err = sum(pred_labels(:) ~= true_labels(:))/length(true_labels);
fprintf('done\n');


fprintf('Training time: %f[s], #trn_errors: %d, nIter: %d\n',...
        ocas_stat.total_time, ocas_stat.nTrnErrors, ocas_stat.nIter);
fprintf('Objval  primal: %f, dual: %f, gap: %f\n', ...
        ocas_stat.Q_P, ocas_stat.Q_D,ocas_stat.Q_P-ocas_stat.Q_D);
fprintf('Testing error: %f %%\n',ocas_tst_err*100);

