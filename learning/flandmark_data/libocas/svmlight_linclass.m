% SVMLIGHT_LINCLASS classifies examples in SVM^light file by linear rule.
%
% Synopsis:
%  [score,true_lab] = svmlight_linclass(data_file,W,[])
%  [score,true_lab] = svmlight_linclass(data_file,W,W0)
%  [score,true_lab] = svmlight_linclass(data_file,W,W0,verb)
% 
% Input:
%  data_file [string] Path to file with examples stored in SVM^light format.
%  W [nDims x nModels] Parameter vectors of nModels linear classifiers.
%  W0 [nModels x 1] Bias of decision rule. If W0 is empty then W0 = 0 is used.
%  verb [1x1] If ~= 0 then prints info (default 0).
%
% Output:
%  score [nModels x nExamples] score(i,j) = W(:,i)'*X_j + W0(i)
%  true_labels [nExamples x 1] labels from the data_file
%
