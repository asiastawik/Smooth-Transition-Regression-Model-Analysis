clear all
close all
clc

load ('Report4_2'); 
K = size(x,2);

LAMBDA = logspace(-2,2,50); %0.01 = 10^(-2), 100 = 10^2, 50 values
C = linspace(quantile(s, 0.1), quantile(s, 0.9), 25); %s - from the data
for i = 1:length(LAMBDA)
    for j = 1:length(C)
        lambda = LAMBDA(i);
        c = C(j);
        G = 1./(1+exp(-lambda*(s-c)));
        Z = x.*G;
        Z2 = x.*(1-G);
        X2 = [Z2 Z];
        betas=OLS(X2,y);
        beta1 = betas(1:K);
        beta2star = betas(K+1:end);
        beta2 = beta2star - beta1;
        params = [beta1; beta2; lambda; c];
        LS(i,j) = loss(x,y,s,params);
    end
end

[ii, jj] = find(LS==min(min(LS)));
lambda = LAMBDA(ii) %23
c = C(jj) %1.5
G = 1./(1+exp(-lambda*(s-c)));
Z = x.*G;
Z2 = x.*(1-G);
X2 = [Z2 Z];
betas=OLS(X2,y);
beta1 = betas(1:K)
beta2star = betas(K+1:end);
beta2 = beta2star - beta1
params = [beta1;beta2;lambda;c];
L = loss(x,y,s,params) %1239.7

disp(['Minimal loss score is in position [ii, jj], where ii = ', num2str(ii), ', jj = ', num2str(jj), ', L is ', num2str(L)])
disp(['Lambda is: ', num2str(lambda)])
disp(['c is: ', num2str(c)])
disp(['beta1 is: ', num2str(beta1')])
disp(['beta2 is: ', num2str(beta2')])