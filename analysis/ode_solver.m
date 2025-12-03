
L = 1200;
v = 600;

k_on = 0.000358;
c = 1.5;
k_off = 0.5; 

syms y(x)
ode = (2*y - 1)*diff(y,x) + k_on*c*(L/v)*(1-y) - k_off*(L/v)*y == 0;
cond1 = y(0) == 0.5;
cond2 = y(1) == 0.5;
conds = [cond1 cond2];

ySol = dsolve(ode, cond1);

fplot(ySol, [0 1])