function [x,stat] = libqp_splx(H,f,b,I,S,x0,opt)
% LIBQP_SPLX Solves quadratic programming task with simplex constraints.
%
% Synopsis:
%  [x,stat] = libqp_splx(H,f,b,I,S)
%  [x,stat] = libqp_splx(H,f,b,I,S,x0)
%  [x,stat] = libqp_splx(H,f,b,I,S,x0,opt)
%  [x,stat] = libqp_splx(H,f,b,I,S,[],opt)
% 
% Description:
%  This function solves the following convex quadratic programing task: 
%
%   min 0.5*x'*H*x + f'*x  
%    x
%
% subject to:   sum(x(find(I==k))) == b(k)  for all k in find(S==0)
%               sum(x(find(I==k))) <= b(k)  for all k in find(S==1)
%                             x(i) >= 0 for all i=1:n
%
% where H is a symmetric positive semi-definite matrix, 
% I is a vector of indices such that unique(I) = 1:max(I)
% and max(I) <= length(I). 
%
% Reference:
%  The algorithm is a generalization of the method proposed in
%   V. Franc, V. Hlavac. A Novel Algorithm for Learning Support Vector Machines
%   with Structured Output Spaces. Research Report K333 22/06, CTU-CMP-2006-04. 
%   May, 2006. ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-TR-2006-04.ps
%  
% Input:
%  H [n x n] Symmetric positive semi-definite matrix. It can be both dense or
%     sparse matrix. The following input vectors, however, must be all dense. 
%  f [n x 1] Vector.
%  b [m x 1] Vector of positive numbers.
%  I [n x 1] Vector of indices such that unique(I) = 1:max(I);
%  S [m x 1] Vector defining constraint sign: 0.. equality, 1..inequality.
%
% Optional inputs: 
%  x0 [n x 1] Initial solution.
%
%  opt [struct] Solver setting:
%   .MaxIter [1 x 1] Maximal number of iterations (default 0xFFFFFFFF);
%   .TolAbs [1 x 1] Absolute tolerance (default 0); halts if QP-QD <= TolAbs  
%   .TolRel [1 x 1] Relative tolerance (default 1e-9); halts if QP-QD <= abs(QP)*TolRel  
%   .QP_TH  [1 x 1] Desired primal value (default -inf); halts if QP <= QP_TH  
%   .verb [1 x 1] If verb > 0 then prints info every verb-th iterations.
%
%  where QP is primal objective value and QD is dual objective value.
%
% Output:
%  x [n x 1] Solution vector.
%
%  stat [struct] Statistics:
%   .QP [1 x 1] Primal objective value.
%   .QD [1 x 1] Dual objective value.
%   .nIter [1 x 1] Number of iterations.
%   .exitflag [1 x 1] Indicates which stopping condition was used:
%     -1  ... Not enough memory.
%      0  ... Maximal number of iterations reached: nIter >= MaxIter.
%      1  ... Relative tolerance reached: QP-QD <= abs(QP)*TolRel
%      2  ... Absolute tolerance reached: QP-QD <= TolAbs
%      3  ... Objective value reached threshold: QP <= QP_TH.
%
% Examples:
%  % setup QP task
%  n=200; X = randn(n,n); H = X'*X; f = randn(n,1); 
%  b = 10*rand(4,1);
%  I=[ones(1,50) 2*ones(1,50) 3*ones(1,50) 4*ones(1,50)]; I=I(randperm(n));
%  S = [0 0 1 1]; 
%
%  % call LIBQP_SPLX solver
%  tic; [x1,stat1] = libqp_splx(H,f,b,I,S); fval1=stat1.QP, toc
%
%  % convert QP task to QUADPROG agruments and run QUADPROG
%  Aequ = zeros(2,n); Aequ(1,find(I==1)) = 1; Aequ(2,find(I==2)) = 1; 
%  bequ = b(1:2);
%  Aneq = zeros(2,n); Aneq(1,find(I==3)) = 1; Aneq(2,find(I==4)) = 1;
%  bneq = b(3:4);
%  tic; [x2,fval2] = quadprog(H,f,Aneq,bneq,Aequ,bequ,zeros(n,1)); fval2, toc
%

% Copyright (C) 2006-2008 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
% Center for Machine Perception, CTU FEL Prague

% options
if nargin < 5, error('At least 5 input arguments are required.'); end

if nargin < 6 || isempty(x0),     

    % create feasible solution
    x0 = zeros(length(f),1); 
    for k=1:length(S)
       if S(k) == 0
         u = min(find(I==k));
         x0(u) = b(k);
       end
    end
end

if nargin < 7, opt = []; end
if ~isfield(opt,'TolAbs'), opt.TolAbs = 0; end
if ~isfield(opt,'TolRel'), opt.TolRel = 1e-9; end
if ~isfield(opt,'MaxIter'), opt.MaxIter = inf; end
if ~isfield(opt,'QP_TH'), opt.QP_TH = -inf; end
if ~isfield(opt,'verb'), opt.verb = 0; end

[x,stat.QP,stat.QD,stat.exitflag,stat.nIter] = ...
    libqp_splx_mex(H,full(diag(H)),f,b,uint32(I),uint8(S),x0,opt.MaxIter,opt.TolAbs,opt.TolRel,opt.QP_TH,opt.verb); 

return;
% EOF