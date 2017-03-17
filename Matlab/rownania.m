function dx = rownania(t,x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)

% x1=x
% x2=x'
% x3=theta
% x4=theta'
% y = [x; x'; o; o']
dx = zeros(4,1);
%c5 = c5*cos(x(3));
%c9 = c9*cos(x(3));
dx(1) = x(2);
dx(3) = x(4);
dx(4) = (c5*cos(x(3))*(c6-c8*x(2)-c10*sin(x(3))*x(4)^2)-c2*c7*x(2)-c3*c7-c4*c7*sin(x(3)))/(c7*c1+c5*c9*cos(x(3))^2);
dx(2) = (c6-c8*x(2)-c9*cos(x(3))*dx(4)-c10*sin(x(3))*x(4)^2)/c7;
%dx(4) = c6/c11 - c7*c2/(c5*c11)*x(2) - c7*c3/(c5*c11) - c7*c4/(c5*c11)*sin(x(3)) - c8/c11*x(2) - c10/c11*sin(x(3))*x(4)^2;
%dx(2) = c1*c6/(c5*c11) - c1*c7*c2/(c5*c5*c11)*x(2) - c1*c7*c3/(c5*c5*c11) - c1*c7*c4/(c5*c5*c11)*sin(x(3)) - c1*c8/(c5*c11)*x(2) - c1*c10/(c5*c11)*sin(x(3))*x(4)^2 + c2/c5*x(2) + c3/c5 + c4/c5*sin(x(3));