% Function to evaluate quadratic basis function & derivative
%   at given quadrature point
%
function[bb,bx,by]=qbf(x,y,it,in,xc,yc,node)
if in<=3,
    in1=in;
	in2=rem(in,3)+1;
	in3=rem(in+1,3)+1;
	i1=node(it,in1);
	i2=node(it,in2);
	i3=node(it,in3);
	d=(xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1));
	t=1+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d;
	bb=t*(2*t-1);
	bx=(yc(i2)-yc(i3))*(4*t-1)/d;
	by=(xc(i3)-xc(i2))*(4*t-1)/d;
  else
    inn=in-3;
	in1=inn;
	in2=rem(inn,3)+1;
	in3=rem(inn+1,3)+1;
	i1=node(it,in1);
	i2=node(it,in2);
	i3=node(it,in3);
	j1=i2;
	j2=i3;
	j3=i1;
	d=(xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1));
	c=(xc(j2)-xc(j1))*(yc(j3)-yc(j1))-(xc(j3)-xc(j1))*(yc(j2)-yc(j1));
	t=1+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d;
	s=1+((yc(j2)-yc(j3))*(x-xc(j1))+(xc(j3)-xc(j2))*(y-yc(j1)))/c;
	bb=4*s*t;
	bx=4*(t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d);
	by=4*(t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d);
end
