function [T,Y] = rk4_back(rhs,x,psiT,tsim,h,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)

Y = zeros(length(tsim/h),length(psiT));
tau = (0:h:tsim)';
dtau=diff(tau);
u=ones(size(dtau));

end

