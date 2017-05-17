function [histU,histQ,typeFlag,i] = bfgs(y0,u0,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho)

    Eps0 = 10e-6;
    maxIt = 200;
    i = 0;
    R = 1;
    
    u = u0;
    histU = zeros(length(u0),maxIt);
    histQ = zeros(1,maxIt);
    
    [Q,dQ] = simDiffEq(y0,u0,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho);
    endFlag = (norm(dQ)<Eps0);
    while(~endFlag)
        if(R == 1)
            W = eye(length(u0));
        else
            r = dQ - g;
            s = u - uPrev;
            W = W + (r*r')/(s'*r) - (W*(s*s')*W)/(s'*W*s);
        end
        for k = 1:16
           if(u(k)==0) dQ(k) = min(0,dQ(k)); 
        end
        d = -W\dQ;
        if(d'*dQ>=0)
            if R==0
                R = 1;
                continue
            else
                typeFlag = 'degenerate'
                endFlag = true;
                break
            end
        end
        g = dQ;
        uPrev = u;
        [u,successFlag] = minimize(d,y0,u,Q,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho,R);
        if(~successFlag)
            if(R == 0)
                R = 1;
                continue
            else
                typeFlag = 'degenerate'
                endFlag = true;
                break
            end
        end
        R = 0;
        [Q,dQ] = simDiffEq(y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho);
        i = i+1;
        histU(:,i) = u;
        histQ(1,i) = Q;
        if(i>=maxIt)
            endFlag = true;
            typeFlag = 'maxIt';
        elseif(norm(dQ)<Eps0)
            endFlag = true;
            typeFlag = 'gradNorm';
        end
    end
end