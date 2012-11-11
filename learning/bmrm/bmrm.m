function [W,stat]= bmrm(data,risk,lambda,options)
% BMRM Bundle Method for regularized Risk Minimization.
%
% Synopsis:
%  [W,stat]= bmrm(data,risk,lambda)
%  [W,stat]= bmrm(data,risk,lambda,options)
%
% Description:
%  Bundle Method for regularized Risk Minimization (BMRM) is an algorithm 
%  for minimization of function
%
%    F(W) = 0.5*lambda*W'*W + R(data,W)
%
%  where R(data,W) is an arbitrary convex function of W. BMRM requires 
%  a function which evaluates R and computes its subgradient at W. 
%  A handle to this function is expected as the second input argument
%  "risk". The calling syntax is
%
%    [Fval,subgrad,data] = risk(data)
%    [Fval,subgrad,data] = risk(data,W)
%
%  which returns function value and subgradient of R(data,W) evaluated 
%  at W. Calling risk(data) without the argument W assumes that 
%  W equals zero vector (this is used to by BMRM to get dimension of
%  parameter vector). The argument data is passed to allow the risk 
%  function to store some auxciliary result if needed. 
%  
% Reference: 
%    Teo et al.: A Scalable Modular Convex Solver for Regularized Risk 
%    Minimization, KDD, 2007
%  
% Inputs:
%  data ... arbitrary data type.
%  risk ... handle of the risk function.
%  lambda [1x1] Regularization parameter.
%
%  options [struct] 
%   .TolRel [1x1] Relative tolerance (default 1e-3). Halt optimization
%      if Fp-Fd <= TolRel*Fp holds.
%   .TolAbs [1x1] Absolute tolerance (default 0). Halt if Fp-Fd <= TolAbs
%      where Fp = F(W) and Fd is a lower bound of F(W_optimal). 
%   .MaxIter [1x1] Maximal number of iterations (default inf). Halt 
%      optimization if nIter >= MaxIter .
%   .verb [1x1] if 1 print progress status (default 0).
%   .BufSize [1x1] Size of cutting plane buffer in MB (default 100).
% 
% Outputs:
%  W [nDim x 1] Solution vector.
%  stat [struct] 
%   .Fp [1x1] Primal objective value.
%   .Fd [1x1] Reduced (dual) objective value.
%   .nIter [1x1] Number of iterations.
%   

% 2010-04-12, Vojtech Franc
    
t0 = cputime;

if nargin < 4, options = []; end
if ~isfield(options,'TolRel'), options.TolRel = 1e-3; end
if ~isfield(options,'TolAbs'), options.TolAbs = 0; end
if ~isfield(options,'MaxIter'), options.MaxIter = inf; end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'BufSize'), options.BufSize = 100; end
if ~isfield(options,'BufSizeCP'), options.BufSizeCP = 0; end
if ~isfield(options,'CleanBuffer'), options.CleanBuffer = inf; end

%-------------------------------------------------------
[R,subgrad,data] = risk(data);

% get paramater space dimension
nDim = length(subgrad);

% computes number of cutting planes from given mega bytes
nCutPlanes = (-nDim*8 + sqrt((nDim*8)^2 +4*8*1024^2 * options.BufSize))/16;
options.BufSize = round(nCutPlanes);

if options.BufSize < options.BufSizeCP 
    options.BufSize = options.BufSizeCP;
end

% inital solution
W = zeros(nDim,1);

if ~issparse(subgrad)
    A = zeros(nDim,options.BufSize);
else
    A = sparse(nDim,options.BufSize);
end
b = zeros(options.BufSize,1);
H = zeros(options.BufSize,options.BufSize);

A(:,1) = subgrad;
b(1) = R;
alpha = [];

nIter = 0;
nCP = 0;
exitflag= -1;

% alloc buffers for meassured statistics
hist_Fd = zeros(options.BufSize+1,1);
hist_Fp = zeros(options.BufSize+1,1);
hist_R = zeros(options.BufSize+1,1);
hist_runtime = zeros(options.BufSize+1,1);

hist_Fd(1) = -inf;
hist_Fp(1) = R+0.5*lambda*norm(W)^2;
hist_R(1) = R;
hist_runtime(1) = cputime-t0;

if options.verb
    fprintf('%4d: tim=%.3f, Fp=%f, Fd=%f, R=%f\n', ...
        nIter, hist_runtime(1), hist_Fp(1), hist_Fd(1), hist_R(1));
end

while exitflag == -1
    nIter = nIter + 1;
    nCP = nCP + 1;
    %    H = A(:,1:nIter)'*A(:,1:nIter)/lambda;
    if nCP > 1,
        H(1:nCP-1,nCP) = full(A(:,1:nCP-1)'*A(:,nCP))/lambda;
        H(nCP,1:nCP-1) = H(1:nCP-1,nCP)';
    end
    H(nCP,nCP) = full(A(:,nCP)'*A(:,nCP))/lambda;
    
    % solve reduced problem
    [alpha,stat] = libqp_splx(H(1:nCP,1:nCP),-b(1:nCP),1,ones(1,nCP),1,[alpha;0]);
    
    W = -A(:,1:nCP)*alpha/lambda;

    nzA = sum(alpha > 0);
        
    [R,subgrad,data] = risk(data,W);
        
    A(:,nCP+1) = subgrad;
    b(nCP+1) = R - A(:,nCP+1)'*W;
    
    Fp = R+0.5*lambda*norm(W)^2;
    Fd = -stat.QP;
                  
    if Fp-Fd<= options.TolRel*abs(Fp)
        exitflag= 1;
    elseif Fp-Fd <= options.TolAbs
        exitflag= 2;    
    elseif nIter >= options.MaxIter
        exitflag= 0;
    end 
    
    hist_runtime(nIter+1) = cputime-t0;
    hist_Fp(nIter+1) = Fp;
    hist_Fd(nIter+1) = Fd;
    hist_R(nIter+1) = R;
    
    if mod(nIter,options.verb) == 0 || exitflag ~= -1
        fprintf('%4d: tim=%.3f, Fp=%f, Fd=%f, (Fp-Fd)=%f, (Fp-Fd)/Fp=%f, R=%f, nCP=%d, nzA=%d\n', ...
                nIter, hist_runtime(nIter+1), Fp, Fd, Fp-Fd,(Fp-Fd)/Fp, R, nCP, nzA);
    end
    
    % clean buffer
    if mod(nIter,options.CleanBuffer) == 0 && nzA < nCP
        old_nCP = nCP;
        idx = find(alpha > 0);
        nCP = length(idx);

        alpha = alpha(idx);
        H(1:nCP,1:nCP) = H(idx,idx);
        A(:,1:nCP) = A(:,idx);        
        b(1:nCP) = b(idx);        
        
        fprintf('done. (old_nCP=%d, new_nCP=%d) \n',old_nCP, nCP);
    end
    
end

stat = [];
stat.Fp = Fp;
stat.Fd = Fd;
stat.nIter = nIter;
stat.hist.Fp = hist_Fp(1:nIter+1);
stat.hist.Fd = hist_Fd(1:nIter+1);
stat.hist.R = hist_R(1:nIter+1);
stat.hist.runtime = hist_runtime(1:nIter+1);

return;
