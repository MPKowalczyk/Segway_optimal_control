function dy = rownania(t,y,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)

% y = [x; x'; o; o']
dy = zeros(4,1);
c5 = c5*cos(y(3));
c9 = c9*cos(y(3));
dy(1) = y(2);
dy(2) = c1*c6/(c5*c11) - c1*c7*c2/(c5*c5*c11)*y(2) - c1*c7*c3/(c5*c5*c11) - c1*c7*c4/(c5*c5*c11)*sin(y(3)) - c1*c8/(c5*c11)*y(2) - c1*c10/(c5*c11)*sin(y(3))*y(4)^2 + c2/c5*y(2) + c3/c5 + c4/c5*sin(y(3));
dy(3) = y(4);
dy(4) = c6/c11 - c7*c2/(c5*c11)*y(2) - c7*c3/(c5*c11) - c7*c4/(c5*c11)*sin(y(3)) - c8/c11*y(2) - c10/c11*sin(y(3))*y(4)^2;