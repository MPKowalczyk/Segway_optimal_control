x1 = 0;
x2 = 0;
x3 = 1;
x4 = 0;
u = @(t)t;
t = 1;
K = 10;
fi_max = 0.52;

%%
H = [ 0, 0,                                                                                                                                                                                                          0,                                                                                                                                                                                                                                                                                                                                        0;
 1, 0,                                                                                                                                                                             -(cos(x3) + 1)/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                              (cos(x3)*(cos(x3) + 1))/(cos(x3)^2 + 1) - 1;
 0, 0, - (cos(x3) + u(t)*sin(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 - x2*sin(x3))/(cos(x3)^2 + 1) - (2*cos(x3)*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2, (cos(x3)*(cos(x3) + u(t)*sin(x3) + x4^2*cos(x3)^2 - x4^2*sin(x3)^2 - x2*sin(x3)))/(cos(x3)^2 + 1) - (sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1) - x4^2*cos(x3) + (2*cos(x3)^2*sin(x3)*(cos(x3)*sin(x3)*x4^2 + sin(x3) + x2*(cos(x3) + 1) - u(t)*(cos(x3) - 1)))/(cos(x3)^2 + 1)^2;
 0, 1,                                                                                                                                                                    -(2*x4*cos(x3)*sin(x3))/(cos(x3)^2 + 1),                                                                                                                                                                                                                                                                                  (2*x4*cos(x3)^2*sin(x3))/(cos(x3)^2 + 1) - 2*x4*sin(x3)];

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