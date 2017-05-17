function [Q,dQdU] = simDiffEq(y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho)
%%Program wykonawczy - jeden krok algorytmu

    %Alokacja pamiêci
    t = zeros(cn(end),1);
    y = zeros(cn(end),size(y0,1));
    psi = zeros(size(y,1),size(y,2));
    dQdU = zeros(size(dtau,1),1);

    y(1,:) = y0;
    t(1) = 0;

    tic
    %Policz x w przód
    for j=1:length(dtau)
        h = dtau(j)/n(j);
        h6 = h/6;
        h2 = h/2;
        for i=cn(j):cn(j+1)-1
            k1 = rhs(t(i),y(i,:)',u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            t_n = t(i)+h2;
            y_n = y(i,:)'+k1*h2;
            k2 = rhs(t_n,y_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            y_n = y(i,:)'+k2*h2;
            k3 = rhs(t_n,y_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            t(i+1) = t(i)+h;
            y_n = y(i,:)'+k3*h;
            k4 = rhs(t(i+1),y_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            y(i+1,:) = (y(i,:)' + (k1 + 2*(k2 + k3) + k4)*h6)'; 
        end
    end
    toc

    %Oblicz wartoœc funkcji celu
    Q = objFun(y(end,:),Vmax,rho);

    psiT = [y(end,4)-1-rho*max(0,y(end,1)-Vmax);0;0;y(end,1);0];

    psi(end,:) = psiT;
    x_middle = zeros(size(y0,1),1);
    t1 = t(end);

    tic
    %Policz psi w ty³
    for j = length(dtau):-1:1
        z = zeros(n(j)+1,1);
        h = dtau(j)/n(j);
        h6 = h/6;
        h2 = h/2;
        for i=cn(j+1):-1:cn(j)+1
            tt = t1;
            t1 = t(i-1);
            d = (y(i,:)'-y(i-1,:)')/h;
            A = [tt^3,tt^2,tt,1;t1^3,t1^2,t1,1;3*tt^2,2*tt,1,0;3*t1^2,2*t1,1,0];
            if det(A) == 0
                x_middle = (y(i,:)'+y(i-1,:)')/2;
            else
                for k=1:size(y0,1)
                    b = [y(i,k)';y(i-1,k)';d(k);d(k)];
                    par = A\b;
                    x_middle(k) = polyval(par',tt-0.5*h);
                end
            end
            [k1,k1z] = rhs_psi(tt,y(i,:)',psi(i,:)',u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            t_n = tt-h2;
            psi_n = psi(i,:)'-k1*h2;
            y_n = x_middle;
            [k2,k2z] = rhs_psi(t_n,y_n,psi_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            psi_n = psi(i,:)'-k2*h2;
            [k3,k3z] = rhs_psi(t_n,y_n,psi_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            psi_n = psi(i,:)'-k3*h;
            [k4,k4z] = rhs_psi(t1,y(i-1,:)',psi_n,u(j),sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2);
            psi(i-1,:) = (psi(i,:)' - (k1 + 2*(k2 + k3) + k4)*h6)';
            z(i-cn(j),1) = z(i+1-cn(j),1) - (k1z + 2*(k2z + k3z) + k4z)*h6;
        end
        dQdU(j) = z(1,1);
    end
    toc
% psi0 = psi(1,:)';
% %SprawdŸ poprawnoœæ psi
% tic
% isGood(Q,psi0,y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho)
% toc
% 
% tic
% isZGood(Q,dQdU,y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho)
% toc
end