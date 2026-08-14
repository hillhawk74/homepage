% M-file to set up geometric information for a rectangular domain
% where we are using triangles
% 

nx2=2*nx-1;
ny2=2*ny-1;
%
% zero arrays
nt=2*(nx-1)*(ny-1);
nn=(nx2-1)*(ny2-2);
node=zeros(nt,6);
indx=zeros(nn,1);
xc=zeros(nn,1);
yc=zeros(nn,1);
wq=zeros(3,1);
xq=zeros(nt,3);
yq=zeros(nt,3);
area=zeros(nt);
%
%  loop over points in y- and x- directions
hx=1/(nx-1);
hy=1/(ny-1);
yy=-hy/2;
xx=-hx/2;
ip=0;
it1=-1;
iu=0; 
for jc=1:ny2
  yy=yy+hy/2;
  xx=-hx/2;
  iy=rem(jc,2);
  for ic=1:nx2
    xx=xx+hx/2;
   	ix=rem(ic,2);
	   ip=ip+1;
   	xc(ip)=xx;
   	yc(ip)=yy;
%
%  set node array
   	if ix==1 & iy==1 & ic~=nx2 & jc~=ny2, 
	     	it1=it1+2;
	     	it2=it1+1;
	     	node(it1,1)=ip;
	     	node(it1,2)=ip+nx2+nx2; 
	     	node(it1,3)=ip+nx2+nx2+2; 
	     	node(it1,4)=ip+nx2 ; 
	     	node(it1,5)=ip+nx2+nx2+1; 
	     	node(it1,6)=ip+nx2+1; 
	     	node(it2,1)=ip;
	     	node(it2,2)=ip+2; 
	   	  node(it2,3)=ip+nx2+nx2+2; 
  	 	  node(it2,4)=ip+1; 
	     	node(it2,5)=ip+nx2+2; 
	     	node(it2,6)=ip+nx2+1;
	   end 
%
%  set indx array to give global node number ; check for type of boundary condition
    if ic ==1 | ic==nx2 
		       indx(ip)=0;
      elseif jc==1 | jc==ny2
         indx(ip)=0;
      else
         iu=iu+1;
         indx(ip)=iu;
    end 
%	  
  end
end 
%
nel=it2;
np=ip;
nunk=iu; 
disp(sprintf( ' no. of elements is  %8i',nel))
disp(sprintf(' no. of nodes is  %8i',np))
disp(sprintf(' no. of unknowns is  %8i',nunk))
%
% set quadrature points for midpoint rule and compute area of triangles
nq=3;
wq(1)=1/3;
wq(2)=wq(1);
wq(3)=wq(1);
for it=1:nel
  ip1=node(it,1);
  ip2=node(it,2);
  ip3=node(it,3);
  x1=xc(ip1);
  x2=xc(ip2);
  x3=xc(ip3);
  y1=yc(ip1);
  y2=yc(ip2);
  y3=yc(ip3);
  xq(it,1)=.5*(x1+x2);
  xq(it,2)=.5*(x2+x3);
  xq(it,3)=.5*(x1+x3);
  yq(it,1)=.5*(y1+y2);
  yq(it,2)=.5*(y2+y3);
  yq(it,3)=.5*(y1+y3);
  area(it)=.5*abs(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2));
end
  
