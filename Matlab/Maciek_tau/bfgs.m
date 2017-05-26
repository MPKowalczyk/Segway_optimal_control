function [tabU, tabQ] = bfgs(iter,eps,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u0)

n = length(u0);
u = u0;
i = 0;
tabU = zeros(iter,n);
tabQ = zeros(iter,n);
W = eye(n);

while (i<iter)
   [Q,dQdU] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
   if (~i) % pomijane w pierwszej iteracji
       s = u - uPrev;
       r = dQdU - g;
       W = W + (r*r')/(s'*r) - (W*(s*s')*W)/(s'*W*s);
   end
   d = -W\dQdU;
   % poszukiwanie na kierunku
   
   i = i+1;
end

end

