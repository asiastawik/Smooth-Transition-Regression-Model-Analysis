clear all
close all
clc

%% Load data
load('Report4_3.mat')
K = size(x,2);

%% Prepare a initial parameters vector with given data
param0 = [beta1;beta2;lam;c];

%% Estimate parameters with ’fminunc’

fun = @(param) loss2(y,x,s,param);
L1 = fun(param0)

[param2, L2] = fminunc(fun, param0);

%% Estimate parameters with ’fmincon’

lb = [-inf(1,2*K), 0.01, quantile(s,0.1)];
ub = [inf(1,2*K), 100, quantile(s, 0.9)];
[param3, L3] = fmincon(fun, param0, [], [], [], [], lb, ub);

%% Estimate parameters with  fminconw with StepTolerance

lb = [-inf(1,2*K), 0.01, quantile(s,0.1)];
ub = [inf(1,2*K), 100, quantile(s, 0.9)];
option = optimoptions("fmincon", 'StepTolerance', 0.05);
[param4, L4] = fmincon(fun, param0, [], [], [], [], lb, ub, [], option);