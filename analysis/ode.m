
L = 10000;
v = 600;

k_on = 0.000358;
c = 1.5;
k_off = 0.5; 

syms y(x)
func = (2*y - 1)*diff(y) + k_on*c*(L/v)*(1-y) -k_off*(L/v)*y == 0;

ySol(x) = dsolve(func);

fplot(ySol, [0 1])