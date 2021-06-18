function f=fun(x,ix3,iX13)
%
global ix3 iX13
% volume fraction of polymer in dilute phase
x3 = ix3;
% volume fraction of nonsolvent in c phase
x1 = x(1);
% volume fraction of solvent in c phase
x2 = 1-x3-x(1);
%
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
%{
% molar volume of DMF [cm^3/mol]
v1 = 82.6;
% molar volume of solvent (THF)
v2 = 79.76;
% molar volume of polymer (PIM-1)
v3 = 53763.44086;
%}
s = v1/v2;
r = v1/v3;
%%% Reference test%%
%s = 0.21;
%r = 0.002;
%%%%
% Interaction parameter nonsolvent(1)-polymer(3)
%
X13 = iX13;
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
% data from THF-PS
g23 = 0.372501868;
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
%
A = -0.0039215691215;
B = 0.06458381077;
C = 0.32598705615;
%}
%
u1 = x1/(x1 + x2);
u2 = x2/(x1 + x2);
g12 = A*u2^2 + B*u2 + C;
dg12du = 2*A*u2 + B;
d2g12du = 2*A; 

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
%{
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
%%%%% Equations for critical point %%%%%
G22 = 1/x1 + s*1/x2 - 2*g12 + 2*(1-2*u2)*dg12du + u1*u2*d2g12du;

G23 = 1/x1 - (g12+X13) + s*g23 + u2*(1-3*u2)*dg12du + u1*u2^2*d2g12du;

G33 = 1/x1 + r*1/x3 - 2*X13 + 2*u2^3*dg12du + u2^3*u1*d2g12du;

G222 = 1/x1^2 - s*1/x2^2 - 6*u2/x2*dg12du+(3-6*u2)*u2/x2*d2g12du;

G223 = 1/x1^2 - 6*u2^2/x2*dg12du + 3*(1-2*u2)*u2^2/x2*d2g12du;

G233 = 1/x1^2 - 6*u2^3/x2*dg12du + (3*u2-6*u2^2)*u2^2/x2*d2g12du;

G333 = 1/x1^2 - r*1/x3^2 + 6*u2^4/x2*dg12du + (3*u2^2 - 2*u2^3)*u2^2/x2*d2g12du;

f = G222*G33^2 - 3*G223*G23*G33 + 3*G233*G23^2 -G22*G23*G333;


