function [uex,uexx,uexy]=exact(x,y)
uex=sin(pi*x)*cos(pi*y);
uexx= pi*cos(pi*x)*cos(pi*y);
uexy=- pi*sin(pi*x)*sin(pi*y);
