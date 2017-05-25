function dydt = rhs(t,y,u,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2)
%Prawa strona równania ró¿niczkoego na x

ratio = u/y(1);
mi = K5*miM*y(3)/((K1+y(3)+y(3)*y(3)/K2)*(K5+y(4)+y(4)*y(4)/K6));
Qp1 = K7*ni1*y(3)/((K3+y(3))*(K7+y(4)));
Qp2 = K8*ni2*y(3)/((K4+y(3))*(K8+y(5)));
dydt = [u; (mi-ratio)*y(2); ratio*(sf-y(3))-(Qp1/Yp1+Qp2/Yp2)*y(2); -ratio*y(4)+Qp1*y(2);-ratio*y(5)+Qp2*y(2)];

end