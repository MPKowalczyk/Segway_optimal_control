%% Model Segway
% Autorzy:
% Marcin Kowalczyk
% Maciej Podsiad³o
% Problem sterowania optymalnego dla pojazdu typu Segway.
clear;
close all;

%% Parametry symulacji
Tsim=10;
fs=1e3;

%% Parametry obiektu sterowania
Mp = 10;
l = 1;
Ip = 10;
km = 0.1;
ke = 0.1;
R = 1;
r = 0.2;
Va = 100;
g = 9.81;
Iw = 0.1;
Mw = 0.5;
x0 = [0;0;0;0];

%% Wspó³czynniki modelu matematycznego
c1 = Mp*l^2 + Ip;
c2 = -2*km*ke/(R*r);
c3 = 2*km/R*Va;  % sterowanie
c4 = Mp*g*l;
c5 = -Mp*l;
c6 = 2*km/(R*r)*Va; % sterowanie
c7 = 2*Mw + 2*Iw/r^2 + Mp;
c8 = 2*km*ke/(R*r^2);
c9 = Mp*l;
c10 = -Mp*l;
c11 = c7*c1/c5 + c9;

%% Symulacja zachowania obiektu
czas = 0:1/fs:Tsim;
[T,Y] = ode45(@rownania,czas,x0,[],c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
[Trk4 Yrk4] = rk4(@rownania,x0,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
%% Wykresy wyniku symulacji
figure(1);
plot(T,Y(:,1));
title('po³o¿enie (t)');
xlabel('czas [s]');
ylabel('po³o¿enie');
grid on;
hold on;
plot(T,Yrk4(:,1));
legend('ode45','rk4');
figure(2);
plot(T,Y(:,2),'r');
title('prêdkoœæ liniowa (t)');
xlabel('czas [s]');
ylabel('prêdkoœæ liniowa');
grid on;
hold on;
plot(T,Yrk4(:,2));
legend('ode45','rk4');
figure(3);
plot(T,Y(:,3));
title('wychylenie (t)');
xlabel('czas [s]');
ylabel('wychylenie');
grid on;
hold on;
plot(T,Yrk4(:,3));
legend('ode45','rk4');
figure(4);
plot(T,Y(:,4),'r');
title('prêdkoœæ k¹towa (t)');
xlabel('czas [s]');
ylabel('prêdkoœæ k¹towa');
grid on;
hold on;
plot(T,Yrk4(:,4));
legend('ode45','rk4');