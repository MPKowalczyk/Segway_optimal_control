function [t,x] = rk4_forw(rhs,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,dtau,n,cn,x0,u)
x=zeros(cn(end),length(x0));
t=zeros(cn(end),1);
for j=1:length(dtau)
	h = dtau(j)/n(j);
    h6 = h/6;
    h3 = h/3;
    h2 = h/2;
    for i=cn(j):cn(j+1)-1
        z=x(i,:);
        k1=rhs(z,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        z1=z+h2*k1;
        k2=rhs(z1,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        z2=z+h2*k2;
        k3=rhs(z2,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        z3=z+h*k3;
        k4=rhs(z3,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u(j));
        x(i+1,:)=z+h3*(k2+k3)+h6*(k1+k4);
        t(i+1)=t(i)+h;
    end
end

end