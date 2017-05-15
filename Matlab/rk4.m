function [T,Y] = rk4(rhs,x0,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)


Y = zeros(length(tsim/h),length(x0));
T = (0:h:tsim)';

hp = h/2;
ht = h/3;
hs = h/6;
for i=1:length(T)
    k1 = rhs(T(i),x0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x1=x0+hp*k1;
    k2 = rhs(T(i)+hp,x1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x2=x0+hp*k2;
    k3 = rhs(T(i)+hp,x2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    x3=x0+h*k3;
    k4 = rhs(T(i)+h,x3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    xnext = x0 + ht*(k2+k3) + hs*(k1+k4);
    Y(i,:) = xnext;
    x0 = xnext; 
end

end