function [u_best] = LineSearchMin(u0,d,stepLen0,wspEksp,wspKontr,maxit,Q0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K)

i=0;
u = u0;
step = stepLen0*d;
uTab = zeros(3,length(u0));
% uTab(1,:) = u0';
[Qnew,grad] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u+step);
if(Qnew<Q0) %ekspansja
    while(i<maxit)
        Qprev = Qnew;
        step = wspEksp*step;
        u = u+step;
        [Qnew,grad] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
        if(Qnew>Qprev)
            break;
        end
    end
else %kontrakcja
        
end


% while(i<maxit)
%     [Qnew,grad] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u+step);
%     
%     if(Qnew<Q0) %ekspansja
%         step = wspEksp*step;
%         [Qnew,grad] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u+step);
%     else %kontrakcja
%         
%     end
%     e
%     i=i+1;
% end

