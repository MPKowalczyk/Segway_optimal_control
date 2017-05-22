function [T,Y] = rk4_dpsi(rhs_a,x,psiT,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)


Y = zeros(length(tsim/h),length(psiT));
T = (0:h:tsim)';

hp = h/2;
ht = h/3;
hs = h/6;
psi = psiT;

for i=length(T):-1:1
    if (i>1)
        xf = x(i-1,:);
    else
        xf = x(i,:);
    end
%     xf = x(i-1,:);
    xp = (x(i,:)+xf)/2;
    
    k1 = rhs_a(psi,x(i,:),T(i),c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    p1 = psi-hp*k1;
    k2 = rhs_a(p1,xp,T(i)-hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    p2 = psi-hp*k2;
    k3 = rhs_a(p2,xp,T(i)-hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    p3 = psi-h*k3;
    k4 = rhs_a(p3,xf,T(i)-h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    psi = psi - ht*(k2+k3) - hs*(k1+k4);
    Y(i,:) = psi;
end

end
