% Problem sterowania optymalnego dla pojazdu typu Segway.
clear all;
close all;
format long e;
format compact;
%% Parametry symulacji
Tsim=6;
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
u_max=20;

%% Wspó³czynniki modelu matematycznego
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
fi_max=pi/18;
K=0.3e3;

%% Symulacja stanu
N = 16;
tau = linspace(0,Tsim,N)';
dtau = diff(tau);
%u = u_max*(2*rand(size(dtau))-1);
u=u_max*ones(size(dtau));
h0 = 0.001;
n = ceil(dtau/h0);
cn = cumsum([1;n]);
x0 = [-5;0;-15*pi/180;0;0];
[t, x] = rk4_tau(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);

%% Równania sprzê¿one
psiT = [-x(end,1:4)';-1]; % pochodna wskaŸnika jakoœci po x(T)
[tpsi, psi, Z] = rk4_dpsi_tau(@dpsi_calc,x,psiT,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
dQdU = Z(cn(1:end-1)+[0;ones(length(cn)-2,1)]);
%% Sprawdzenie równañ sprzê¿onych
de=1e-7;
epsilon=de*eye(5);
dQ=zeros(5,1);
Q0=cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
for i=1:5
    dQ(i)=(cost(@rownania_penalty,x0+epsilon(:,i),dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)-Q0)/de;
end
dQ
psi(1,:)'
roznica_sprzezone = dQ + psi(1,:)'
%% Sprawdzenie gradientów
eps = 1e-7;
Q0=cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
dQdU_check = zeros(size(u));
for i=1:length(u)
    du = zeros(size(u));
    du(i)=eps;
    dQdU_check(i) = (cost(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u+du) - Q0)/eps;
end
roznica_gradienty = dQdU - dQdU_check

%% Optymalizacja
iter=200;
e0=1e-8;
tic;
[u,tabU,tabQ]=bfgs2(iter,e0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u);
toc;
%% wykres wychylenia po zastosowaniu sterowania opt
[t, x] = rk4_tau(@rownania_penalty,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
figure(1);
subplot(3,2,1);
plot(t,x(:,1));
xlabel('Czas [t]');
ylabel('Po³o¿enie [m]');
grid on;
% figure(2);
subplot(3,2,2);
plot(t,x(:,2));
xlabel('czas [t]');
ylabel('prêdkoœæ liniowa [m/s]');
grid on;
% figure(3);
subplot(3,2,3);
plot(t,x(:,3));
hold on;
plot(t,fi_max*ones(size(t)),'r',t,-fi_max*ones(size(t)),'r');
xlabel('Czas [t]');
ylabel('Wychylenie [rad]');
grid on;
% figure(4);
subplot(3,2,4);
plot(t,x(:,4));
xlabel('Czas [t]');
ylabel('Prêdkoœæ k¹towa [rad/s]');
grid on;
% figure(5);
subplot(3,2,5);
plot(t,x(:,5));
xlabel('Czas [t]');
ylabel('Kara za wychylenie');
grid on;
% figure(6);
subplot(3,2,6);
for i=1:N-1
    ster(cn(i):cn(i+1)) = ones(length(cn(i):cn(i+1)),1)*u(i);
end
plot(t,ster);
xlabel('Czas [t]');
ylabel('Sterowanie [V]');
grid on;