warning('off','all')

%Ustaw wartoœci parametrów
Vmax = 5;
sf = 330;
miM = 0.679;
K1 = 7.24;
K2 = 351;
K3 = 9.3;
K4 = 62.1;
K5 = 47.7;
K6 = 10.9;
K7 = 171;
K8 = 230;
ni1 = 1.18;
ni2 = 0.141;
Yp1 = 0.492;
Yp2 = 0.213;

% rho = 1000;
rho = 0;

% u = (0:0.025:0.375)';
u = zeros(16,1);
y0 = [1.5; 1.5; 8.6; 0; 0];

step = 0.001;
tau = (0:1:16)';
dtau = diff(tau);
n = ceil(dtau/step);
cn = cumsum([1;n]);

tic
[histU,histQ,typeFlag,itNo] = bfgs(y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho);
toc

histQ(1,1:itNo)