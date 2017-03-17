%% Model Segway
% Autorzy:
% Marcin Kowalczyk
% Maciej Podsiad³o
% Problem sterowania optymalnego dla pojazdu typu Segway.
clear;
close all;

%% Parametry symulacji
Tsim=10;
fs=1e5;

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
%[T,Y] = ode45(@ss_sympy,czas,x0,[],Mp,Mw,l,Ip,Iw,km,ke,R,r,g,Va);

%% Wykresy wyniku symulacji
figure(1);
plot(T,Y(:,1));
title('po³o¿enie (t)');
xlabel('czas [s]');
ylabel('po³o¿enie');
grid on;
figure(2);
plot(T,Y(:,2),'r');
title('prêdkoœæ liniowa (t)');
xlabel('czas [s]');
ylabel('prêdkoœæ liniowa');
grid on;
figure(3);
plot(T,180/pi*Y(:,3));
title('wychylenie (t)');
xlabel('czas [s]');
ylabel('wychylenie');
grid on;
figure(4);
plot(T,Y(:,4),'r');
title('prêdkoœæ k¹towa (t)');
xlabel('czas [s]');
ylabel('prêdkoœæ k¹towa');
grid on;