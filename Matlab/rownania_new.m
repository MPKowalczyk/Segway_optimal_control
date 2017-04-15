function dx = rownania_new(t,x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,u)

a1 = c5*c6*cos(x(3))/c7 - c3;
a2 = -(c5*c8*cos(x(3))/c7 + c2);
a3 = -c5*c10*cos(x(3))*sin(x(3))/c7;
a4 = -c4*sin(x(3));
L = c1 + c5*c9*cos(x(3))^2/c7;

dx = zeros(4,1);
dx(1) = x(2);
dx(3) = x(4);
dx(4) = (a1*u(t) + a2*x(2) + a3*x(4)^2 + a4)/L;
dx(2) = (c6*u(t) - c8*x(2) - c9*cos(x(3))*dx(4) - c10*x(4)^2*sin(x(3)))/c7;
