x1 = 0;
x2 = 0;
x3 = 1;
x4 = 0;
u = @(t)t;
t = 1;
K = 10;
fi_max = 0.52;

%%
% H = [ 0, 0,                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                        0;
%  1, 0,                                                                                                                                                                             -(cos(x3) + 1)/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                              (cos(x3)*(cos(x3) + 1))/(cos(x3)^2 + 1) - 1;
%  0, 0, - (cos(x3) + u(t)*sin(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 - x2*sin(x3))/(cos(x3)^2 + 1) - (2*cos(x3)*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2, (cos(x3)*(cos(x3) + u(t)*sin(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 - x2*sin(x3)))/(cos(x3)^2 + 1) - (sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1) - x4^2*cos(x3) + (2*cos(x3)^2*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2;
%  0, 1,                                                                                                                                                                    -(2*x4*cos(x3)*sin(x3))/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                  (2*x4*cos(x3)^2*sin(x3))/(cos(x3)^2 + 1) - 2*x4*sin(x3)];

 
H =[ 0, 0,                                                                                                                                                                                                                                                                 0,                                                                                                                                                                                                                        0;
 1, 0,                                                                                                                                                                                                                     -(81*cos(x3) + 648/25)/(4*(5*cos(x3)^2 - 16)),                                                                                                                                                                              (81*(cos(x3) + 10))/(20*(5*cos(x3)^2 - 16));
 0, 0, (7848*cos(x3) + 500*x4^2*(2*cos(x3)^2 - 1) - 450*u(t)*sin(x3) + 2025*x2*sin(x3))/(500*cos(x3)^2 - 1600) + (8*cos(x3)*sin(x3)*((25*cos(x3)*sin(x3)*x4^2)/4 + (981*sin(x3))/10 - x2*((405*cos(x3))/16 + 81/10) + u(t)*((45*cos(x3))/8 + 9/5)))/(5*cos(x3)^2 - 16)^2, ((26487*cos(2*x3))/10 + 390*x4^2*cos(x3) + 810*x2*sin(2*x3) + (81*x2*sin(3*x3))/4 + 50*x4^2*cos(3*x3) - 180*sin(2*x3)*u(t) - (9*sin(3*x3)*u(t))/2 + sin(x3)*((5589*x2)/20 - (621*u(t))/10) - 981/2)/(5*cos(2*x3) - 27)^2;
 0, 1,                                                                                                                                                                                                                              -(5*x4*sin(2*x3))/(5*sin(x3)^2 + 11),                                                                                                                                                                                       (20*x4*sin(x3))/(5*sin(x3)^2 + 11)];
 
H1 = psi_value([x1;x2;x3;x4;0],t,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
%% uwzglêdnienie funkcji f5 karz¹c¹ za wychylenia segwaya wiêksze od fi max funkcj¹ kwadratow¹ ze wspó³czynnikiem K
df5 = zeros(4,1);
if( abs(x3)>fi_max )
    df5(3) = K*(abs(x3) - fi_max);
else
    df5(3) = 0;
end

H(:,5) = df5;
H(5,:) = zeros(1,5);


% Z = [ 0, 0,                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                          0;
%  1, 0,                                                                                                                                                                                         -(cos(x3) + 1)/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                                                (cos(x3)*(cos(x3) + 1))/(cos(x3)^2 + 1) - 1;
%  0, 0, - (cos(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 + conj(u)*sin(x3) - x2*sin(x3))/(cos(x3)^2 + 1) - (2*cos(x3)*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - conj(u)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2, (cos(x3)*(cos(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 + conj(u)*sin(x3) - x2*sin(x3)))/(cos(x3)^2 + 1) - (sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - conj(u)*(cos(x3) - 1)))/(cos(x3)^2 + 1) - x4^2*cos(x3) + (2*cos(x3)^2*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - conj(u)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2;
%  0, 1,                                                                                                                                                                                -(2*x4*cos(x3)*sin(x3))/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                                    (2*x4*cos(x3)^2*sin(x3))/(cos(x3)^2 + 1) - 2*x4*sin(x3)]