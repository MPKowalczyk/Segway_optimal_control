function dx = rownania_ss(t,x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,u)

dx = zeros(4,1);
dx(1) = x(2);
dx(2) = -c8/c7*x(2)-c9/c7*((c5*c6/c7*cos(x(3))-c3)*u(t)-(c5*c8/c7*cos(x(3))+c2)*x(2)-c4*sin(x(3))-c5*c10/c7*x(4)^2*sin(x(3))*cos(x(3)))*cos(x(3))/(c1+c5*c9/c7*cos(x(3))^2)-c10/c7*x(4)^2*sin(x(3))+c6/c7*u(t);
dx(3) = x(4);
dx(4) = ((c5*c6/c7*cos(x(3))-c3)*u(t)-(c5*c8/c7*cos(x(3))+c2)*x(2)-c4*sin(x(3))-c5*c10/c7*x(4)^2*sin(x(3))*cos(x(3)))/(c1+c5*c9/c7*cos(x(3))^2);