%function to determine value of u at the boundary given the (x,y) coordinates of the
%point on the boundary
%
function[ub]=ubdry(x,y)
ub=sin(pi*x)*cos(pi*y);
