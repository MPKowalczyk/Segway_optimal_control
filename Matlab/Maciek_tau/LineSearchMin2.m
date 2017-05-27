function [u_best,success] = LineSearchMin2(u0,d,stepLen0,wspEksp,wspKontr,maxit,Q0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A)

i=1;
u = u0;
step = stepLen0*d;
success=0;
%uTab = zeros(3,length(u0));
% uTab(1,:) = u0';
[Qnew,~] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u+step);
if(Qnew<Q0) %ekspansja
    while(i<maxit)
        %disp(['Ekspansja - iteracja: ' num2str(i)]);
        Qprev = Qnew;
        step = wspEksp*step;
        u = u+step;
        [Qnew,~] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u);
        if(Qnew>Qprev)
            success=1;
            break;
        end
        i=i+1;
    end
else %kontrakcja
    while(i<maxit)
        %disp(['Kontrakcja - iteracja: ' num2str(i)]);
        step=wspKontr*step;
        u=u+step;
        [Qnew,~] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,A,u);
        if(Qnew<Q0)
            success=1;
            break;
        end
        i=i+1;
    end
end
if success
    u_best=u;
else
    u_best=u0;
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