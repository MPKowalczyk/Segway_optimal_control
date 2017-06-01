function [dpsi, dz] = dpsi_calc(psi,x,t,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u)

[H, df_u]=psi_value(x,t,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,fi_max,K,u);
dpsi = -H*psi;
dz = psi'*df_u;

end

