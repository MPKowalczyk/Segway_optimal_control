function H = psi_value(x,t,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)
k1=-c8/c7;
k2=c9/c7;
k3=c10/c7;
k4=c6/c7;
k5=-c5*c6/c7;
k6=-c3;
k7=c5*c8/c7;
k8=-c2;
k9=c4;
k10=-c5*c10/c7;
k11=c1;
k12=c5*c9/c7;

L=(k5*cos(x(3))+k6)*u(t)+(k7*cos(x(3))+k8)*x(2)+k9*sin(x(3))+k10*x(4)^2*sin(x(3))*cos(x(3));
M=k11+k12*cos(x(3))^2;
S=L/M;
dL_x3=k10*(2*cos(x(3))^2-1)*x(4)^2+k9*cos(x(3))-k5*u(t)*sin(x(3)) - k7*x(2)*sin(x(3));
dL_x4=k10*x(4)*sin(2*x(3));
dM_x3=-k12*sin(2*x(3));

dS_x3=dL_x3/M-L/M^2*dM_x3;
dS_x4=dL_x4/M;
if abs(x(3))>fi_max
    df5_x3=K*(abs(x(3))-fi_max);
else
    df5_x3=0;
end

H=...
[ 0,                                                         0, 0,     0,      0;...
  1,                                                        k1, 0,     0,      0;...
  0, k3*x(4)^2*cos(x(3)) + k2*cos(x(3))*dS_x3 - k2*sin(x(3))*S, 0, dS_x3, df5_x3;...
  0,                  2*k3*x(4)*sin(x(3)) + k2*cos(x(3))*dS_x4, 1, dS_x4,      0;...
  0,                                                         0, 0,     0,      0];
end