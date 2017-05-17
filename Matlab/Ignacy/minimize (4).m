function [uOpt,successFlag] = minimize(d,y0,u,Q,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho,R)
    
    if(R)
        EpsQ = 10e-15;
    else
        EpsQ = 10e-6;
    end

    uOpt = u + d;
    [newQ,~] = simDiffEq(y0,uOpt,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho);
    if(newQ<Q)
        successFlag = true;
        return
    end
    while(newQ>=Q && norm(uOpt-u)>EpsQ)
        d = 0.5*d;
        uOpt = u + d;
        [newQ,~] = simDiffEq(y0,uOpt,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho);
    end
    if(newQ>=Q)
        uOpt = u;
        successFlag = false;
    else
        successFlag = true;
    end
end