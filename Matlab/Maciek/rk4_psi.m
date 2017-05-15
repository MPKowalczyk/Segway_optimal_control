function [T,Y] = rk4_psi(rhs,x,psiT,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)


Y = zeros(length(tsim/h),length(psiT));
T = (0:h:tsim)';
h=-h;
% T = (tsim:h:0)';

hp = h/2;
ht = h/3;
hs = h/6;
psi = psiT;

% for i=1:length(T)
for i=length(T):-1:1
%     odwr = length(T)+1-i;
    odwr = i;
    k1 = rhs(T(i),x(odwr,:),psi,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
%     x1=x0+hp*k1;
    p1 = psi+hp*k1;
%     k2 = rhs(T(i)+hp,x1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,u);
    k2 = rhs(T(i)+hp,x(odwr,:),p1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
%     x2=x0+hp*k2;
    p2 = psi+hp*k2;
%     k3 = rhs(T(i)+hp,x2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,u);
    k3 = rhs(T(i)+hp,x(odwr,:),p2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
%     x3=x0+h*k3;
    p3 = psi+h*k3;
%     k4 = rhs(T(i)+h,x3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,u);
    k4 = rhs(T(i)+h,x(odwr,:),p3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    psi = psi + ht*(k2+k3) + hs*(k1+k4);
    Y(odwr,:) = psi; % pytanie co z pierwsz� iteracj�?
%     x0 = xnext; 
end

% T = flip(T);

end