from sympy import collect, Eq, pprint, solve, symbols

# Counting

x1,x2,x3,x4,x5=symbols('x1 x2 x3 x4 x5')
Mp,Mw,l,Ip,Iw,km,ke,R,r,Va,g,sin,cos=symbols('Mp Mw l Ip Iw km ke R r Va g sin cos')
eq1=Eq((Mp*l**2+Ip)*x5-2*km*ke*x1/R/r+2*km*Va/R+Mp*g*l*sin,-Mp*l*x2*cos)
eq2=Eq(2*km*Va/R/r,(2*Mw+2*Iw/r**2+Mp)*x2+2*km*ke*x1/R/r**2+Mp*l*x5*cos-Mp*l*x4**2*sin)
sol=solve([eq1,eq2],x2,x5)
sol_x2=collect(sol.get(x2),[x1,x2,x3,x4,x5])
sol_x5=collect(sol.get(x5),[x1,x2,x3,x4,x5])

# Printing
pprint(sol)