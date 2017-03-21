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
km = 0.9;
ke = 0.9;
R = 1;
r = 0.2;
Va = 0;
g = 9.81;
Iw = 0.1;
Mw = 0.5;
x0 = [0;0;pi+10*pi/180;0];

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
t = 0:1/fs:Tsim;
[t,x] = ode45(@rownania,t,x0,[],c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);

%% Wykresy wyniku symulacji
figure(1);
plot(t,x(:,1));
title('Po³o¿enie');
xlabel('Czas [s]');
ylabel('x_1 [m]');
grid on;
figure(2);
plot(t,x(:,2));
title('Prêdkoœæ liniowa');
xlabel('Czas [s]');
ylabel('x_2 [m/s]');
grid on;
figure(3);
plot(t,(x(:,3)-pi)*180/pi);
title('Wychylenie');
xlabel('Czas [s]');
ylabel('x_3 [^o]');
grid on;
figure(4);
plot(t,x(:,4)*180/pi);
title('Prêdkoœæ k¹towa');
xlabel('Czas [s]');
ylabel('x_4 [^o/s]');
grid on;