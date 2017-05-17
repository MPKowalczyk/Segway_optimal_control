function [dpsidt,dz] = rhs_psi(t,y,psi,u,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2)
%Prawa strona równania ró¿niczkoego na psi

ratio = u/y(1);
ratio2 = ratio/y(1);
mi = K5*miM*y(3)/((K1+y(3)+y(3)*y(3)/K2)*(K5+y(4)+y(4)*y(4)/K6));
Qp1 = K7*ni1*y(3)/((K3+y(3))*(K7+y(4)));
Qp2 = K8*ni2*y(3)/((K4+y(3))*(K8+y(5)));

a12 = -y(2)*ratio2;
a13 = ratio2*(sf-y(3));
a14 = -y(4)*ratio2;
a15 = -y(5)*ratio2;
a22 = ratio - mi;
a23 = Qp1/Yp1 + Qp2/Yp2;
a24 = -Qp1;
a25 = -Qp2;
a32 = -y(2)*mi*(K1 - (y(3)^2/K2))/(y(3)*(K1+y(3)+(y(3)^2/K2)));
a34 = -y(2)*K3*Qp1/(y(3)*(K3+y(3)));
a35 = -y(2)*K4*Qp2/(y(3)*(K4+y(3)));
a33 = (ratio - a34/Yp1 - a35/Yp2);
a42 = y(2)*mi*(1+2*y(4)/K6)/(K5+y(4)+(y(4)^2/K6));
a43 = -y(2)*Qp1/(Yp1*(K7+y(4)));
a44 = ratio + y(2)*Qp1/(K7+y(4));
a53 = -y(2)*Qp2/(Yp2*(K8+y(5)));
a55 = ratio + y(2)*Qp2/(K8+y(5));

A = [0,a12,a13,a14,a15;0,a22,a23,a24,a25;0,a32,a33,a34,a35;0,a42,a43,a44,0;0,0,a53,0,a55];
dpsidt = A*psi;
dz = psi(1)-psi(2)*y(2)/y(1)+psi(3)*(sf-y(3))/y(1)-psi(4)*y(4)/y(1)-psi(5)*y(5)/y(1);

end