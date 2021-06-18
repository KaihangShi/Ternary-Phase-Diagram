clear
clc
global ix3 iX13

options = optimoptions('fmincon','MaxIterations',100000,'OptimalityTolerance',1e-10);
% Might be good choice for initial guess
% Initial guess for DMF-THF-PIM1 system
%x0 = [0.2;0.5;0.5];
% Initial guess for Toluene-THF-PIM1 system
x0 = [0.6];
%lb = [0];
%ub = [1];

%A = [1];
%b = [1];

% Input parameters
nloop = 50;
%nloop = 1;
% polymer volume fraction in dilute phase (no need to change except want to 
% cover larger range)
ix3 = 1e-5;
% interaction parameter nonsolvent(1)-polymer(3)
iX13 =1.7;
%xsol = zeros(2*nloop,4);
xsol = zeros(nloop,4);
base = 1e-5;

for i = 1:nloop
    %
    if mod(i,10) == 0
        base = base*10;
    end
    ix3 = ix3+base;
    %}
    %ix3d = ix3d*0.1;
    %[x,f] = fmincon(@fun,x0,A,b,[],[],lb,ub,[],options);
    %j = i*2;
    %xsol(j-1,1) = fun(x);
    %xsol(j-1,2) = x(1);
    %xsol(j-1,4) = x(2);
    %xsol(j-1,3) = 1-x(1)-x(2);
    [x,fval] = fsolve(@fun,x0);
    xsol(i,1) = fval;
    xsol(i,2) = x(1);
    xsol(i,4) = 1- x(1)-ix3;
    xsol(i,3) = ix3;
end
