function [u,tabU, tabQ] = bfgs2(iter,eps,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u0)

n = length(u0);
u = u0;
i = 1;
tabU = zeros(iter,n);
tabQ = zeros(iter,1);
R=1;
stepLen0=1e-0;
wspEksp=4;
wspKontr=0.9;
[Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u);
dQdU(u==u_max & dQdU<0)=0;
dQdU(u==-u_max & dQdU>0)=0;
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
%    disp(['Iloczyn skalarny: ' num2str(-d'*dQdU)]);
   
%    [u,success]=LineSearchMin2(u,d,stepLen0,wspEksp,wspKontr,maxit,Q,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K);
   [u_best,Q_best,EkspSuccess] = ekspansja(u,d,stepLen0,wspEksp,wspKontr,maxit,Q,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max);
   u = u_best(2,:)';
%    if(~EkspSuccess)
%        [u_best,Q_best,KontrSuccess] = kontrakcja(u,d,stepLen0,wspEksp,wspKontr,maxit,Q,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K);
%        disp('kontrakcja');
%    end
%    success = EkspSuccess || KontrSuccess;
   
   success = EkspSuccess;
%    u = u_best(2,:)';
    u(u>u_max)=u_max;
    u(u<-u_max)=-u_max;
   
   g=dQdU;
   Qprev = Q;
   [Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u);
   dQdU(u==u_max & dQdU<0)=0;
   dQdU(u==-u_max & dQdU>0)=0;
   if(Q>Qprev)
       disp('B³¹d metody, pogorszenie funkcji kosztu! Iteracja:');
       disp(i);
       break;
   end
   tabU(i+1,:)=u;
   tabQ(i+1)=Q;
   if success
       R=0;
   elseif ~R
       R=1;
   else
       tabU=tabU(1:i+1,:);
       tabQ=tabQ(1:i+1);
       disp('Brak poprawy na kierunku najszybszego spadku');
       break;
   end
   i = i+1;
%    disp(['Nowy koszt: ' num2str(Q) ' | Norma gradientu: ' num2str(norm(dQdU))]);
disp([Q norm(dQdU)]);
end

end