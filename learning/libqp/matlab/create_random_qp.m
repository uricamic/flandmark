% This script generates randomly several testing 
% Quadratic Programming tasks.
%
% The generated tasks are stored to ../data/ folder.
%

randn('state',0);
rand('state',0);

%=======================================================
% Quadratic Programming task with simplex constraints
% solved by LIBQP_SPLX
%=======================================================

% TASK 1: splxcon_qp_random100x10
% number of variables 100
% number of equlity and inequality constraints 10 (not considering constraints x >= 0)
I = [];
for i=1:10;
  I = [I i*ones(1,10)];
end
Z = randn(length(I),length(I)*2); H=Z*Z';
f=randn(length(I),1); 
b = rand(10,1)*10;

% 5 randomly selected constraints will be equlities and inequlities
S = zeros(10,1);
idx = randperm(10);
S(idx(1:5)) = 1;

save ../data/splx_qp_random100x10 H f b I S;

% TASK 2: splxcon_qp_random300x30
% number of variables 300
% number of equlity and inequality constraints 30 (not considering constraints x >= 0)
I = [];
for i=1:30;
  I = [I i*ones(1,10)];
end
Z = randn(length(I),length(I)*2); H=Z*Z';
f=randn(length(I),1); 
b = rand(30,1)*10;

% 5 randomly selected constraints will be equlity and inequlity
S = zeros(30,1);
idx = randperm(30);
S(idx(1:15)) = 1;

save ../data/splx_qp_random300x30 H f b I S;

%=======================================================
% Quadratic Programming task solved by GSMO algorithm
% solved by LIBQP_GSMO
%=======================================================

% TASK 1: gsmo_qp_random100
% number of variables 100

f=randn(100,1)*20; 
Z = randn(length(f),length(f)*2); H=Z*Z';

a = randn(1,length(f)); a(find(a==0)) = 1;
b = randn(1);
UB = rand(length(f),1);
LB = zeros(length(f),1);

save ../data/gsmo_qp_random100 H f a b UB LB;

% TASK 2: gsmo_qp_random300
% number of variables 300

f=randn(300,1)*20; 
Z = randn(length(f),length(f)*2); H=Z*Z';

a = randn(1,length(f)); a(find(a==0)) = 1;
b = randn(1);
UB = zeros(length(f),1);
LB = -rand(length(f),1);

save ../data/gsmo_qp_random300 H f a b UB LB;
