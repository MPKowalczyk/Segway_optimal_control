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
x0 = [0;0;10*pi/180;0;0];

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

%% Sprawdzenie
de=1e-7;
epsilon=de*eye(5);
dQ=zeros(5,1);
Q0=cost(@rownania_penalty,x0,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,@(t) 0);
for i=1:5
    dQ(i)=(cost(@rownania_penalty,x0+epsilon(:,i),Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,@(t) 0)-Q0)/de;
end

%% Wyliczenie Psi
[trk,xrk] = rk4(@rownania_penalty,x0,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,@(s) 0);
figure(1);
plot(trk,xrk(:,3));
psiT = [-xrk(end,1:4)';-1]; % pochodna wska�nika jako�ci po x(T) DO ZMODYFIKOWANIA bez f5
[tpsi, psi] = rk4_dpsi(@dpsi_calc,xrk,psiT,Tsim,1/fs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,@(s) 0);
dQ
psi0 = psi(1,:)
roznica = dQ + psi(1,:)'