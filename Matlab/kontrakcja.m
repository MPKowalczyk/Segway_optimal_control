function [u_best,Q_best,success] = kontrakcja(u0,d,stepLen0,wspEksp,wspKontr,maxit,Q0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max)


Qlastlast = Q0;
Qlast = Q0;
step = stepLen0*d;
uPrevPrev = u0;
uPrev = u0;
success=0;
u_best = u0;
Q_best = Q0;

for i=1:maxit
    u = u0+step;
    %     [Qnew,grad] = testowa(u);
    [Qnew,~] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u_max,u);
    if(Qnew-Q0<0)
        disp('Poprawa kontrakcja');
        success = 1;
        u_best = u;
        Q_best = Qnew;
        break;
    end
    step = wspKontr*step;
    uPrevPrev = uPrev;
    uPrev = u;
    Qlastlast = Qlast;
    Qlast=Qnew;
end
% u_best = [u uPrev uPrevPrev]';
% Q_best = [Qnew Qlast Qlastlast]';
end


