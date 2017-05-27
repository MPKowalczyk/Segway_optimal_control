function [u_best,Q_best,success] = ekspansja(u0,d,stepLen0,wspEksp,wspKontr,maxit,Q0,x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K)

Qlastlast = Q0;
Qlast = Q0;
step = stepLen0*d;
uPrevPrev = u0;
uPrev = u0;
success = 0;

for i=1:maxit
    u = u0+step;
%     [Qnew,~] = testowa(u);
    [Qnew,~] = SimSegway(x0,dtau,cn,h0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
    if(Qnew-Qlast>1e-8)
        break;
    end
    step = wspEksp*step;
    uPrevPrev = uPrev;
    uPrev = u;
    Qlastlast = Qlast;
    Qlast=Qnew;
    success=1;
end

u_best = [uPrevPrev uPrev u]';
Q_best = [Qlastlast Qlast Qnew]';
end
