clear all
close all
clc

t = linspace(-2*pi(),2*pi(),500)';
s = cos(t);
lambda = 3;
c = 0.3;

G = 1-exp(-lambda*(s-c).^2);

X = [cos(linspace(-1/2,1/2*pi(),500)') ones(500,1)]; %1/4*2pi = 1/2pi
beta1 = [1 ; 2];
beta2 = [2 ; 1];

y = X*beta1 + X.*G*beta2 + randn(500,1);
figure(1)
plot(t, y)
hold on
plot(t, G)
hold off
title("Smooth transition regression and transition function.")
legend('y', 'G')

loss_value = loss(X,y,s,[beta1;beta2;lambda;c])