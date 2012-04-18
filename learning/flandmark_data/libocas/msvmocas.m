% MSVMOCAS Train multi-class linear SVM classifier using OCAS solver.
%
% Synopsis:
%  [W,stat] = msvmocas(X,y,C,Method,TolRel,TolAbs,QPBound,BufSize,nExamples,MaxTime,verb)
%  [W,stat] = msvmocas(svmlight_data_file,C,Method,TolRel,TolAbs,QPBound,BufSize,nExamples,MaxTime,verb)
%
% Desription:
%  This function trains multi-class linear SVM classifier by solving
%
%      W^* = argmin 0.5*sum_y (W(:,y)'*W(:,y)) + C*  sum     max( (y~=y(i)) + (W(:,y) - W(:,y(i))'*X(:,i))
%              W                                   i=1:nData   y
%
%  The function accepts examples either in Matlab matrix X (both sparse and dense) and 
%  a dense vector y or as path to a file in SVM^light format.
%
% Reference:
%  V. Franc, S. Sonnenburg. To be published. 2009
%
% Input:
%   data_file [string] Training examples stored in SVM^light format.
%
%   X [nDim x nExamples] training inputs (sparse or dense matrix).
%   y [nExamples x 1] labels; intgers 1,2,...nY
%   C [1x1] regularization constant
%   Method [1x1] 0..cutting plane; 1..OCAS  (default 1)
%   TolRel [1x1] halts if Q_P-Q_D <= abs(Q_P)*TolRel  (default 0.01)
%   TolAbs [1x1] halts if Q_P-Q_D <= TolAbs  (default 0)
%   QPValue [1x1] halts if Q_P <= QPBpound  (default 0)
%   BufSize [1x1] Initial size of active constrains buffer (default 2000)
%   nExamples [1x1] Number of training examplesused for training; must be >0 and <= size(X,2).
%     If nExamples = inf then nExamples is set to size(X,2).
%   MaxTime [1x1] halts if time used by solver (data loading time is not counted) exceeds
%    MaxTime given in seconds. Use MaxTime=inf (default) to switch off this stopping condition. 
%   verb [1x1] if non-zero then prints some info; (default 1)
%
% Output:
%   W [nDim x nY] Paramater vectors of decision rule; [dummy,ypred] = max(W'*x)
%   stat [struct] Optimizer statistics (field names are self-explaining).
%
% Example:
%  C = 1; 
%
%  % case 1: loading data directly from file in SVM^light format
%  [W,stat] = msvmocas('example4_train.light',C);
%
%  % case 2: using data loaded in Matlab
%  load('example4_train','X','y');
%  [W,stat] = msvmocas(X,y,C);
%
%  % classification
%  load('example4_test.mat','X','y');
%  [dummy,ypred] = max(W'*X);
%  sum(ypred(:) ~= y(:))/length(y)
% 

%
% Copyright (C) 2008 Vojtech Franc, xfrancv@cmp.felk.cvutr.cz
%                    Soeren Sonnenburg, soeren.sonnenburg@first.fraunhofer.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public 
% License as published by the Free Software Foundation; 
% Version 3, 29 June 2007
