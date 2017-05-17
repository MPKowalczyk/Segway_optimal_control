function Q = objFun( xT, Vmax, rho )
%Funkcja celu

Q = xT(1)*(1-xT(4)) + 0.5*rho*max(0,xT(1)-Vmax);

end