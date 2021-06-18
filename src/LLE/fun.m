function f=fun(x)
%

% volume fraction of polymer in dilute phase
x3d = 1e-20;
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
%
% molar volume of nonsolvent [cm^3/mol]
v1 = 106.3;
% molar volume of solvent (THF)
v2 = 79.76;
% molar volume of polymer (PIM-1)
v3 = 53763.44086;
% RT [J/mol]
RT = 2454.01281;
s = v1/v2;
r = v1/v3;
% Interaction parameter nonsolvent-PS
X13 = 1.6;
% Interaction parameter THF-PS
X23 = 0.6;

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

mu1c = log(x1c) + 1 - x1c - s*x2c - r*x3c +...
    (g12c*x2c + X13*x3c)*(x2c+x3c) - ...
    s*X23*x2c*x3c - x2c*u1c*u2c*dgdu2c;
mu1d = log(x1d) + 1 - x1d - s*x2d - r*x3d +...
    (g12d*x2d + X13*x3d)*(x2d+x3d) - ...
    s*X23*x2d*x3d - x2d*u1d*u2d*dgdu2d;
mu2c = s*log(x2c) + s - x1c - s*x2c - r*x3c...
    + (g12c*x1c + X23*s*x3c)*(x1c+x3c) - ...
    X13*x1c*x3c + x1c*u1c*u2c*dgdu2c;
mu2d = s*log(x2d) + s - x1d - s*x2d - r*x3d...
    + (g12d*x1d + X23*s*x3d)*(x1d+x3d)...
    - X13*x1d*x3d + x1d*u1d*u2d*dgdu2d;
mu3c = r*log(x3c) + r - x1c - s*x2c - r*x3c...
     + (X13*x1c + s*X23*x2c)*(x1c + x2c) -... 
    g12c*x1c*x2c;
mu3d = r*log(x3d) + r - x1d - s*x2d - r*x3d...
     + (X13*x1d + s*X23*x2d)*(x1d + x2d) - g12d*x1d*x2d;

f = (mu1c - mu1d)^2 + (mu2c - mu2d)^2 + (mu3c - mu3d)^2;


