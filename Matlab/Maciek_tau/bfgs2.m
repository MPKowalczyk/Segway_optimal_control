function [tabU, tabQ] = bfgs2(iter,eps,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u0)

n = length(u0);
u = u0;
i = 1;
tabU = zeros(iter,n);
tabQ = zeros(iter,1);
R=1;
stepLen0=1;
wspEksp=1.1;
wspKontr=0.9;
[Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
tabQ(1)=Q;
tabU(1,:)=u;

while (i<iter && norm(dQdU)>eps)
   if(R)
       W=eye(n);
       maxit=100;
   else
       maxit=10;
       s = u - uPrev;
       r = dQdU - g;
       W = W + (r*r')/(s'*r) - (W*(s*s')*W)/(s'*W*s);
   end
   d = -W\dQdU;
   uPrev=u;
   [u,success]=LineSearchMin2(u,d,stepLen0,wspEksp,wspKontr,maxit,Q,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K);
   g=dQdU;
   [Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
   tabU(i+1,:)=u;
   tabQ(i+1)=Q;
   if success
       R=0;
   elseif ~R
       R=1;
   else
       tabU=tabU(1:i+1,:);
       tabQ=tabQ(1:i+1);
       break;
   end
   i = i+1;
end

end