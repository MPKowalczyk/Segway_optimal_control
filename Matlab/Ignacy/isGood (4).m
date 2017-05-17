function isPsiGood = isGood(Q,psi0,y0,u,Vmax,dtau,n,cn,sf,miM,K1,K2,K3,K4,K5,K6,K7,K8,ni1,ni2,Yp1,Yp2,rho)
%Sprawdzenie poprawnoœci psi

Eps = 10e-8;

%Alokacja pamiêci
t = zeros(cn(end),1);
y = zeros(cn(end),size(y0,1));
dQ = zeros(size(y0,1),1);

for k=1:size(y0)
    x_delta = y0;
    x_delta(k) = x_delta(k) + Eps;
    y(1,:) = x_delta;
    
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
    
    deltaQk = objFun(y(end,:),Vmax,rho);
    dQ(k) = (deltaQk - Q)/Eps;
end

% [dQ,-psi0]
isPsiGood = ((norm(dQ + psi0)/norm(psi0))<0.005);

end