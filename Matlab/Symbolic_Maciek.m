syms x1 x2 x3 x4 a1(x3) a2(x3) a3(x3) a4(x3) u(t) L(x3) f1(x1,x2,x3,x4) f2(x1,x2,x3,x4) f3(x1,x2,x3,x4) f4(x1,x2,x3,x4) real
syms psi1 psi2 psi3 psi4

c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
c5 = 1;
c6 = 1;
c7 = 1;
c8 = 1;
c9 = 1;
c10 = 1;

a1 = c5*c6*cos(x3)/c7 - c3;
a2 = -(c5*c8*cos(x3)/c7 + c2);
a3 = -c5*c10*cos(x3)*sin(x3)/c7;
a4 = -c4*sin(x3);
L = c1 + c5*c9*cos(x3)^2/c7;

f1 = x2;
f2 = x4;
f3 = (a1*u(t) + a2*x2 + a3*x4^2 + a4)/L;
f4 = (c6*u(t) - c8*x2 - c9*cos(x3)*f3 - c10*x4^2*sin(x3))/c7;

J = jacobian([f1; f2; f3; f4], [x1 x2 x3 x4]);
H = J.';
dpsi=-H*[psi1; psi2; psi3; psi4];