function [R,subgrad,data] = risk_svm(data,W)
% RISK_SVM evaluates SVM risk term given by L1-hinge loss.
% 
% Synopsis:
%  [R,subgrad,data] = risk_svm(data)
%  [R,subgrad,data] = risk_svm(data,W)
%
% Description:
%  Let the risk term be defined as
%
%  R(W) = 1/m sum_{i=1}^m max(0, 1-data.y(i)*(W'*[data.X(:,i);data.X0]))
%
%  This function returns value R and subgradient SUBGRAD of the 
%  risk R(W) at W.
%

if nargin < 2
    % evaluates the risk and the subgradient of at W = 0

    nExamples = size(data.X,2);
    R = 1;
    if data.X0 == 0
        subgrad = -data.X*data.y(:)/nExamples;
    else        
        subgrad = [-data.X*data.y(:)/nExamples; ...
                   -data.X0*sum(data.y)/nExamples];
    end
               
else
        
    nExamples = size(data.X,2);
        
    if data.X0 == 0
        dfce = 1 - (W'*data.X).*reshape(data.y,1,nExamples);

        idx = find(dfce>0);                
        R = sum(dfce(idx))/nExamples;
        
        subgrad = -data.X(:,idx)*reshape(data.y(idx),length(idx),1)/nExamples;
    else
    
        dfce = 1 - (W(1:end-1)'*data.X + W(end)*data.X0).*(data.y(:)');
        idx = find(dfce>0);
        
        R = sum(dfce(idx))/nExamples;
        subgrad = [-data.X(:,idx)*reshape(data.y(idx),length(idx),1)/nExamples; ...
                   -data.X0*sum(data.y(idx))/nExamples];
    end
        
end
   
return;
% EOF