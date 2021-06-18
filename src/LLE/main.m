clear
clc

options = optimoptions('fmincon','MaxIterations',100000,'OptimalityTolerance',1e-10);
% Might be good choice for initial guess
%x0 = [0.1;0.5;0.499];
x0 = [0.6;0.3;0.6]
lb = [0;0;0];
ub = [1;1;1];

A = [0 0 0;0 0 0;1 1 0];
b = [0;0;1];

[x,f] = fmincon(@fun,x0,A,b,[],[],lb,ub,[],options)

fun(x)