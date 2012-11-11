function [ R, subgrad, data ] = risk2( data, W )
%RISK evaluates risk term.
% 
% Synopsis:
%  [R,subgrad,data] = risk(data)
%  [R,subgrad,data] = risk(data,W)
%
% Description:
%  Let the risk term be defined as
%
%  R(W) = 1/m sum_{i=1}^m [ max_{y \in Y}( (L(y_i, y) + <w, \Psi(x_i, y)> ) - <w, \Psi(x_i, y+i)> ]
%
%  This function returns value R and subgradient SUBGRAD of the 
%  risk R(W) at W.
%    
% 10-08-10 Michal Uricar
% 11-07-11 Michal Uricar, corners dataset (checked-only) + loss function modification

    options = data.options;

    if (nargin < 2)
        W = buildPsi2(options, data.tmpData, data.Y{1});
        W = sparse(length(W), 1);
    end
    
    psi_xiyhat = sparse(length(W), data.nImages);
    psi_xiyi = sparse(length(W), data.nImages);
    suma_all = zeros(1, data.nImages);
    
    for i = 1 : data.nImages
%     parfor i = 1 : data.nImages
        [lbpdat lbp_sparse] = getPsiMat(data, i);
        GT = data.Y{i};
        
        % \hat{y_i} = argmax_{y \in Y} [ L(y_i, y) + <w, \Psi(x_i, y)> ]
        L = computeL(options, GT, data.kappa(i));       % with normalization
%         L = computeL(options, GT);                      % without normalization
        y_hat = argmax_mex(options, W, data.mapTable, lbp_sparse, L);
        
        % \Psi(x_i, \hat{y_i})
        psi_xiyhat(:, i) = buildPsi2(options, lbpdat, y_hat);
        % \Psi(x_i, y_i)
        psi_xiyi(:, i) = buildPsi2(options, lbpdat, GT);
        
        % sum_{i = 1}^{m} ( max_{y \in Y}( L(y_i, y) + <w, \Psi(x_i, y_i)> )> ) - <w, \Psi(x_i, y_i)> )
        % with normalization coefficient kappa        
        suma_all(i) = 100 * data.kappa(i) * ...
                        sum(1/data.options.M * sqrt(sum( (y_hat - GT).^2 ))) ...
                        + W'*psi_xiyhat(:, i) - W'*psi_xiyi(:, i);
%         % without normalization coefficient kappa
%         suma_all(i) = sum(1/data.options.M * sqrt(sum( (y_hat - GT).^2 ))) ...
%                         + W'*psi_xiyhat(:, i) - W'*psi_xiyi(:, i);
    end;

    suma = sum(suma_all);
    
    % \hat{R}(w) = 1/m \sum_{i = 1}^m [ max_{y \in Y}( L(y_i, y) + <w, \Psi(x_i, y_i)> )> ) - <w, \Psi(x_i, y_i)> ]
    R = 1/data.nImages * suma;
    % \partial_w \hat{R}(w) = 1/m \sum_{i = 1}^m ( \Psi(x_i, \hat{y_i}) - \Psi(x_i, y_i) )
    subgrad = 1/data.nImages * sum(psi_xiyhat - psi_xiyi, 2);

end
