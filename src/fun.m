function f=fun(x,ix3d,iX13)
%
global ix3d iX13
% volume fraction of polymer in dilute phase
x3d = ix3d;
% volume fraction of nonsolvent in c phase
x1c = x(1);
% volume fraction of solvent in c phase
x2c = x(2);
%
x3c = 1-x1c-x2c;
% volume fraction of nonsolvent in dilute phase
x1d = x(3);
% volume fraction of solvent in dilute phase
x2d = 1-x1d-x3d;
%{
% molar volume of nonsolvent [cm^3/mol]
% Toluene
v1 = 106.3;
% molar volume of solvent
% THF
v2 = 79.76;
% molar volume of polymer
% PIM1
v3 = 100000;
%v3 = 53763.44086;
%}
%
% molar volume of DMF [cm^3/mol]
v1 = 82.6;
% molar volume of solvent (THF)
v2 = 79.76;
% molar volume of polymer (PIM-1)
v3 = 100000;
%}
s = v1/v2;
r = v1/v3;
%%% Reference test%%
%s = 0.21;
%r = 0.002;
%%%%
% Interaction parameter nonsolvent(1)-polymer(3)
%
X13c = iX13;
X13d = iX13;
dX13c = 0;
dX13d = 0;
%}
%{
k1 = 0.568;
k2 = -0.5583;
k3 = 0.4722;
X13c = k1*x3c^2 + k2*x3c + k3;
X13d = k1*x3d^2 + k2*x3d + k3;
dX13c = 2*k1*x3c + k2;
dX13d = 2*k1*x3d + k2;
%}
%%% Reference test %%
%%X13 = 1;
%%%
% Interaction parameter solvent(2)-polymer(3)
%
kk = 0;
% data from THF-PS
%bb = 0.4;
bb = 0.372501868;
%bb = 0.789773848;
g23c = kk*x3c + bb;
g23d = kk*x3d + bb;
dg23c = kk;
dg23d = kk;
%}
%%%% Parameter derived from vapor pressure
%{
kk = -0.4828;
bb = 0.2261;
g23c = kk*x3c + bb;
g23d = kk*x3d + bb;
dg23c = kk;
dg23d = kk;
%}
%%% data from literature
%{
a1 = 1.396;
b1 = -2.1027;
c1 = 0.573;
g23c = a1*x3c^2 + b1*x3c + c1;
g23d = a1*x3d^2 + b1*x3d + c1;
dg23c = 2*a1*x3c + b1;
dg23d = 2*a1*x3d + b1;
%}
%{
% Interaction parameter nonsolvent(1)-solvent(2)
A = 1.297;
B = -0.9707;
C0 = -0.06669;
%}
%{
A = -0.3782;
B = 0.3712;
C0 = 0.1708;
%}
%{
A = 0;
B = 0;
C0 = 0.1708;
%}
%{
u1c = x1c/(x1c + x2c);
u2c = x2c/(x1c + x2c);
u1d = x1d/(x1d + x2d);
u2d = x2d/(x1d + x2d);
dgdu2c = B*C0/(1 - C0*u2c)^2;
dgdu2d = B*C0/(1 - C0*u2d)^2;
g12c = A + B/(1 - C0*u2c);
g12d = A + B/(1 - C0*u2d);
%}
% Polynimial fit VLE
%
%{
A = -0.0039215691215;
B = 0.06458381077;
C = 0.32598705615;
u1c = x1c/(x1c + x2c);
u2c = x2c/(x1c + x2c);
u1d = x1d/(x1d + x2d);
u2d = x2d/(x1d + x2d);
dgdu2c = 2*A*u2c + B;
dgdu2d = 2*A*u2d + B;
g12c = A*u2c^2 + B*u2c + C;
g12d = A*u2d^2 + B*u2d + C;
%}
% LLE
%{
p1 = -0.05437;
p2 = 0.1953;
p3 = -0.3768;
p4 = 0.5063;
p5 = -1.318;
u1c = x1c/(x1c + x2c);
u2c = x2c/(x1c + x2c);
u1d = x1d/(x1d + x2d);
u2d = x2d/(x1d + x2d);
g12c = p1*u2c^4 + p2*u2c^3 + p3*u2c^2 + p4*u2c + p5;
g12d = p1*u2d^4 + p2*u2d^3 + p3*u2d^2 + p4*u2d + p5;
dgdu2c = 4*p1*u2c^3 + 3*p2*u2c^2 + 2*p3*u2c + p4;
dgdu2d = 4*p1*u2d^3 + 3*p2*u2d^2 + 2*p3*u2d + p4;
%}
%
% Interaction parameter nonsolvent (DMF) - solvent (THF)
% parameters for g12 fitting polynomial
p1 = 0.5777;
p2 = -0.5866;
p3 = 0.4738;
p4 = 0.2363;
p5 = 0.5821;
% Parameters using normal Q (THF) for UNIFAC
%p1 = 1.028;
%p2 = -1.149;
%p3 = 0.8065;
%p4 = 0.2597;
%p5 = 0.7186;
u1c = x1c/(x1c + x2c);
u2c = x2c/(x1c + x2c);
u1d = x1d/(x1d + x2d);
u2d = x2d/(x1d + x2d);
g12c = p1*u2c^4 + p2*u2c^3 + p3*u2c^2 + p4*u2c + p5;
g12d = p1*u2d^4 + p2*u2d^3 + p3*u2d^2 + p4*u2d + p5;
dgdu2c = 4*p1*u2c^3 + 3*p2*u2c^2 + 2*p3*u2c + p4;
dgdu2d = 4*p1*u2d^3 + 3*p2*u2d^2 + 2*p3*u2d + p4;
%}
%%%%% Equations for phase equilibira %%%%%
mu1c = log(x1c) + 1 - x1c - s*x2c - r*x3c +...
    (g12c*x2c + X13c*x3c)*(x2c+x3c) - ...
    s*g23c*x2c*x3c - x2c*u1c*u2c*dgdu2c - ...
    x1c*x3c^2*dX13c - s*x2c*x3c^2*dg23c;

mu1d = log(x1d) + 1 - x1d - s*x2d - r*x3d +...
    (g12d*x2d + X13d*x3d)*(x2d+x3d) - ...
    s*g23d*x2d*x3d - x2d*u1d*u2d*dgdu2d - ...
    x1d*x3d^2*dX13d - s*x2d*x3d^2*dg23d;

mu2c = s*log(x2c) + s - x1c - s*x2c - r*x3c...
    + (g12c*x1c + g23c*s*x3c)*(x1c+x3c) - ...
    X13c*x1c*x3c + x1c*u1c*u2c*dgdu2c - ...
    x1c*x3c^2*dX13c - s*x2c*x3c^2*dg23c;

mu2d = s*log(x2d) + s - x1d - s*x2d - r*x3d...
    + (g12d*x1d + g23d*s*x3d)*(x1d+x3d) - ...
    X13d*x1d*x3d + x1d*u1d*u2d*dgdu2d - ...
    x1d*x3d^2*dX13d - s*x2d*x3d^2*dg23d;

mu3c = r*log(x3c) + r - x1c - s*x2c - r*x3c...
     + (X13c*x1c + s*g23c*x2c)*(x1c + x2c) -... 
     g12c*x1c*x2c + (x1c*dX13c + s*x2c*dg23c)*x3c*(x1c+x2c);

mu3d = r*log(x3d) + r - x1d - s*x2d - r*x3d...
     + (X13d*x1d + s*g23d*x2d)*(x1d + x2d) - ...
     g12d*x1d*x2d + (x1d*dX13d + s*x2d*dg23d)*x3d*(x1d+x2d);

f = (mu1c - mu1d)^2 + (mu2c - mu2d)^2 + (mu3c - mu3d)^2;


