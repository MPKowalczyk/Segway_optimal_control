function [T,Y] = rk4_tau(rhs,x0,dtau,cn,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)

Y = zeros(cn(end),length(x0));
T = zeros(cn(end),1);

hp = h/2;
ht = h/3;
hs = h/6;

x = x0;
Y(1,:)=x;
for j=1:length(dtau)
    for i=cn(j):cn(j+1)-1
        k1 = rhs(T(i),x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        x1=x+hp*k1;
        k2 = rhs(T(i)+hp,x1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        x2=x+hp*k2;
        k3 = rhs(T(i)+hp,x2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        x3=x+h*k3;
        k4 = rhs(T(i)+h,x3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        x = x + ht*(k2+k3) + hs*(k1+k4);
        Y(i+1,:) = x;
        T(i+1) = T(i)+h;
    end
end

end