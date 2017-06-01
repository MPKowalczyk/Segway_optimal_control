function [T,Y,Z] = rk4_dpsi_tau(rhs_a,x,psiT,dtau,cn,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)


Y = zeros(size(x));
T = zeros(cn(end),1);
Z = zeros(cn(end),1);

hp = h/2;
ht = h/3;
hs = h/6;
ho = h/8;
psi = psiT;
Y(end,:)=psi;

for j=length(dtau):-1:1
    dxp = rownania_penalty(0,x(cn(j+1),:),c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
    Z(cn(j+1))=0;
    for i=cn(j+1):-1:cn(j)+1    
        [k1,zk1] = rhs_a(psi,x(i,:),T(i),c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        p1 = psi-hp*k1;
        dx1 = rownania_penalty(0,x(i-1,:),c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        xp = (x(i-1,:) + x(i,:))'/2 + (dx1-dxp)*ho;
        dxp = dx1;
        [k2,zk2] = rhs_a(p1,xp,T(i)-hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        p2 = psi-hp*k2;
        [k3,zk3] = rhs_a(p2,xp,T(i)-hp,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        p3 = psi-h*k3;
        [k4,zk4] = rhs_a(p3,x(i-1,:),T(i)-h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        psi = psi - ht*(k2+k3) - hs*(k1+k4);
        Z(i-1) = Z(i) - ht*(zk2+zk3) - hs*(zk1+zk4);
        Y(i-1,:) = psi;
%         Z(i)   = z;
    end
    dQdU(j) = Z(cn(j)+1);
end

end
