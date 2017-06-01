function [Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u)
%% Ograniczenie sterowania
u(u>u_max)=u_max;
u(u<-u_max)=-u_max;
[t, x] = rk4_tau(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
Q = 0.5*(x(end,1:4)*x(end,1:4)')+x(end,5);
%% Równania sprzê¿one
psiT = [-x(end,1:4)';-1]; % pochodna wskaŸnika jakoœci po x(T)
[tpsi, psi, Z] = rk4_dpsi_tau(@dpsi_calc,x,psiT,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
dQdU = Z(cn(1:end-1)+[0;ones(length(cn)-2,1)]);

end