function [T Y] = rk4(rhs,x0,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)


Y = zeros(length(tsim/h),length(x0));
T = (0:h:tsim)';

hp = h/2;
ht = h/3;
hs = h/6;
for i=1:length(T)
    k1 = rhs(T(i),x0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
    k2 = rhs(T(i)+hp*k1,x0+hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
    k3 = rhs(T(i)+hp*k2,x0+hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
    k4 = rhs(T(i)+h*k3,x0+h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10);
    xnext = x0 + ht*(k2+k3) + hs*(k1+k4);
    Y(i,:) = xnext;
    x0 = xnext; 
end

end