%% Calculate derivatives for optimal control
syms x1 x2 x3 x4 x5 real
syms k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 real
syms u real
L=(k5*cos(x3)+k6)*u+(k7*cos(x3)+k8)*x2+k9*sin(x3)+k10*x4^2*sin(x3)*cos(x3);
M=k11+k12*cos(x3)^2;
syms S(x3,x4) f5(x3,x4)
f1=x2;
f2=k1*x2+k2*S*cos(x3)+k3*x4^2*sin(x3)+k4*u;
f3=x4;
f4=S;

%% Calculate derivatives
H=jacobian([f1,f2,f3,f4,f5],[x1,x2,x3,x4,x5]).';