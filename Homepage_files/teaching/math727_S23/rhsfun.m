% function for the right-hand side of problem
%
function[rhs]=rhsfun(x,y)
rhs=2*pi*pi*sin(pi*x)*cos(pi*y);
