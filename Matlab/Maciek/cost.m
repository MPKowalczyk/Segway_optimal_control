function Q = cost(rhs,x,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)

T = (0:h:tsim)';

hp = h/2;
ht = h/3;
hs = h/6;
for i=1:length(T)
    k1 = rhs(T(i),x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x1=x+hp*k1;
    k2 = rhs(T(i)+hp,x1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x2=x+hp*k2;
    k3 = rhs(T(i)+hp,x2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x3=x+h*k3;
    k4 = rhs(T(i)+h,x3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x = x + ht*(k2+k3) + hs*(k1+k4);
end
% Wyliczeni kosztu na podstawie x(T)
Q=0.5*(x(1:4)'*x(1:4))+x(5);
%%
%x(3) niestabilny punkt pionowy wedlug tej notacji to 3.14! Powinno byæ 0 !

end