function dx = rownania_penalty(t,x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)

dx = zeros(5,1);
dx(1) = x(2);
dx(3) = x(4);
dx(4) = (-c5*cos(x(3))*(c6*u(t)-c8*x(2)+c10*sin(x(3))*x(4)^2)-c2*c7*x(2)-c3*u(t)*c7+c4*c7*sin(x(3)))/(c7*c1+c5*c9*cos(x(3))^2);
dx(2) = (c6*u(t)-c8*x(2)+c9*cos(x(3))*dx(4)+c10*sin(x(3))*x(4)^2)/c7;
if abs(x(3))-fi_max>0
    dx(5)=K*(abs(x(3))-fi_max)^2/2;
else
    dx(5)=0;
end