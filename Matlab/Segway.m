%% Model Segway
% Autorzy:
% Marcin Kowalczyk
% Maciej Podsiad³o
% Problem sterowania optymalnego dla pojazdu typu Segway.

%% Parametry symulacji
Tsim=10;
fs=1e4;

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
fi_max=pi/6;
K=1e1;

%% Sterowanie optymalne
N=11;
tau=linspace(0,Tsim,N)';
dtau = diff(tau);
n = ceil(dtau*fs);
cn = cumsum([1;n]);
u=ones(size(dtau));
x0 = [0;0;10*pi/180;0;0];

%% Symulacja zachowania obiektu
%t = 0:1/fs:Tsim;
%tic;
[t,x] = rk4_forw(@rownania_penalty,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,dtau,n,cn,x0',u);
%T1=toc;
%tic;
[trk,xrk] = rk4(@rownania_penalty,x0,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,1);
%[t_ode,x_ode] = ode45(@rownania_penalty,t,x0,[],c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,1);
%T2=toc;
%Q=cost(@rownania_penalty,x0,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,@(t) 0);

%% Wykresy wyniku symulacji
figure(1);
plot(t,x(:,1),'b',trk,xrk(:,1),'r');
title('Po³o¿enie');
xlabel('Czas [s]');
ylabel('x_1 [m]');
grid on;
legend('New','Old');
figure(2);
plot(t,x(:,2),'b',trk,xrk(:,2),'r');
title('Prêdkoœæ liniowa');
xlabel('Czas [s]');
ylabel('x_2 [m/s]');
grid on;
legend('New','Old');
figure(3);
plot(t,x(:,3),'b',trk,xrk(:,3),'r');
title('Wychylenie');
xlabel('Czas [s]');
ylabel('x_3 [rad]');
grid on;
legend('New','Old');
figure(4);
plot(t,x(:,4),'b',trk,xrk(:,4),'r');
title('Prêdkoœæ k¹towa');
xlabel('Czas [s]');
ylabel('x_4 [rad/s]');
grid on;
legend('New','Old');