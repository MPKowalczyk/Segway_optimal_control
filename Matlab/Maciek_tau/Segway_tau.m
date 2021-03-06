% Problem sterowania optymalnego dla pojazdu typu Segway.
clear all;
close all;

%% Parametry symulacji
Tsim=1;
fs=1e3;

%% Parametry obiektu sterowania
Mp = 10;
l = 1;
Ip = 10;
km = 0.9;
ke = 0.9;
R = 1;
r = 0.2;
Va = 0;
g = 9.81;
Iw = 0.1;
Mw = 0.5;

%% Wsp�czynniki modelu matematycznego
c1 = Mp*l^2 + Ip;
c2 = -2*km*ke/(R*r);
c3 = 2*km/R;%*Va;  % sterowanie
c4 = Mp*g*l;
c5 = -Mp*l;
c6 = 2*km/(R*r);%*Va; % sterowanie
c7 = 2*Mw + 2*Iw/r^2 + Mp;
c8 = 2*km*ke/(R*r^2);
c9 = Mp*l;
c10 = -Mp*l;
c11 = c7*c1/c5 + c9;
fi_max=pi/6;
K=1e3;

%% Symulacja stanu
A=eye(4);
%A(3,3)=10;
N = 5;
tau = linspace(0,Tsim,N)';
dtau = diff(tau);
u = ones(size(dtau));
h0 = 0.001;
n = ceil(dtau/h0);
cn = cumsum([1;n]);
x0 = [0;0;10*pi/180;0;0];
[t, x] = rk4_tau(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);

%% R�wnania sprz�one
%psiT = [-x(end,1:4)';-1]; % pochodna wska�nika jako�ci po x(T)
psiT = [-A*x(end,1:4)';-1];
[tpsi, psi, Z] = rk4_dpsi_tau(@dpsi_calc,x,psiT,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
dQdU = Z(cn(1:end-1)+[0;ones(length(cn)-2,1)]);
%% Sprawdzenie r�wna� sprz�onych
de=1e-7;
epsilon=de*eye(5);
dQ=zeros(5,1);
Q0=cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u);
for i=1:5
    dQ(i)=(cost(@rownania_penalty,x0+epsilon(:,i),dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u)-Q0)/de;
end
dQ
psi(1,:)'
roznica_sprzezone = dQ + psi(1,:)'
%% Sprawdzenie gradient�w
eps = 1e-7;
Q0=cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u);
dQdU_check = zeros(size(u));
for i=1:length(u)
    du = zeros(size(u));
    du(i)=eps;
    dQdU_check(i) = (cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u+du) - Q0)/eps;
end
roznica_gradienty = dQdU - dQdU_check

%% Optymalizacja
iter=100;
e0=1e-8;
[tabU,tabQ]=bfgs2(iter,e0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u);
[t, x] = rk4_tau(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,tabU(end,:));
figure(1);
plot(tabQ);
title('Optimization process');
xlabel('Iteration');
ylabel('Cost function');
grid on;
figure(2);
plot(t,x(:,1));
title('Po�o�enie');
xlabel('Czas [s]');
ylabel('x_1 [m]');
grid on;
figure(3);
plot(t,x(:,2));
title('Pr�dko��');
xlabel('Czas [s]');
ylabel('x_2 [m/s]');
grid on;
figure(4);
plot(t,x(:,3));
title('Nachylenie');
xlabel('Czas [s]');
ylabel('x_3 [rad]');
grid on;
figure(5);
plot(t,x(:,4));
title('Pr�dko�� k�towa');
xlabel('Czas [s]');
ylabel('x_4 [rad/s]');
grid on;