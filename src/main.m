clear
clc
global ix3d iX13

options = optimoptions('fmincon','MaxIterations',100000,'OptimalityTolerance',1e-10);
% Might be good choice for initial guess
% Initial guess for DMF-THF-PIM1 system
%x0 = [0.2;0.5;0.5];
% Initial guess for Toluene-THF-PIM1 system
x0 = [0.1;0.7;0.25];
lb = [0;0;0];
ub = [1;1;1];

A = [0 0 0;0 0 0;1 1 0];
b = [0;0;1];

% Input parameters
nloop = 1;
% polymer volume fraction in dilute phase (no need to change except want to 
% cover larger range)
ix3d = 1e-70;
% interaction parameter nonsolvent(1)-polymer(3)
iX13 =1.7;
xsol = zeros(2*nloop,4);
%base = 1e-20;

for i = 1:nloop
    %{
    if mod(i,10) == 0
        base = base*10;
    end
    ix3d = ix3d+base;
    %}
    %ix3d = ix3d*0.1;
    [x,f] = fmincon(@fun,x0,A,b,[],[],lb,ub,[],options);
    j = i*2;
    xsol(j-1,1) = fun(x);
    xsol(j-1,2) = x(1);
    xsol(j-1,4) = x(2);
    xsol(j-1,3) = 1-x(1)-x(2);
    xsol(j,1) = fun(x);
    xsol(j,2) = x(3);
    xsol(j,4) = 1- x(3)-ix3d;
    xsol(j,3) = ix3d;
end
