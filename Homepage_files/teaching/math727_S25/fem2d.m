% c *****************************************************************
%c
%c  This program solves -(laplacian u) =f in a box in R2 where 
%c  Dirichlet   boundary conditions are imposed.      The
%c  code uses continuous piecewise quadratic basis functions on
%c  triangles.
%c
%c *****************************************************************
%c *****************************************************************
%c
%c  INPUT
%c    nx    -  number of points in the x-direction
%c    ny    -  number of points in the y-direction
%c    xl,xr -  left and right x-coordinates of domain
%c    yb,yt -  bottom and top y-coordinates of domain
%c    nlocn-  number of local nodes per element (for quadratics it's 6)
%c    nq    -  number of quadrature points per element used in the assembly
%c    nqe   -  number of quadrature points per element used in the error routine
%c    lda   -  leading dimension of a(i,j) 
%c    lde   -  leading dimension of arrays dimensioned by number of elements 
%c
%c *****************************************************************
%c
%c  GEOMETRY
%c    xc(i)  -  x-coordinate of global node number i
%c    yc(i)  -  y-coordinate of global node number i
%c    node(it,j)  - array giving global node number corresponding to
%c                  local node number j of element number it
%c    area(it)    - area of element it
%c    indx(i)     - gives unknown number at global node number i
%c                  (i=0 if no unknown at node)
%c    xq(it,i)    - x-coordinate of quadrature point j for element it
%c    yq(it,i)    - y-coordinate of quadrature point j for element it
%c
%c    nunk        - number of unknowns
%c    nel         - number of elements
%c    np          - number of global nodes (or points)
%c *****************************************************************
%c
%c  ASSEMBLY
%c    a(i,j)  - coefficient matrix stored here using banded storage
%c              dimensioned by the number of unknowns by total bandwidth
%c    f(i)    - right hand side array which is overwritten by solver with
%c              solution
%c
%c *****************************************************************
%c *****************************************************************
%c
%c  USER SUPPLIED ROUTINES
%c    exact - gives the exact solution of your PDE (used for error calculation)
%c    rhs - gives the right hand side of your differential equation  
%c    ubdry - gives the value of the your solution on the boundary
%c
%c *****************************************************************
%c
% NOTE - for efficiency we should be using a banded Cholesky instead of the built-in
% matlab function "slash"
% *****************************************************************
% 
% calculate solution for two grids
disp(' solution of -laplacian u = f,  u=0 on boundary of (0,1)x(0,1)')
nlocn=6;
for ngrid=1:2
   if (ngrid==1)
        nx=5;
     elseif(ngrid==2)
        nx=9;
  end 
  ny=nx;
  disp(' ')
  disp(sprintf('number of points in x- and y-directions is %8i',nx))
%
% call geometry routine to set up geometry arrays and
%   quadrature information
  geom
%
%  assemble coefficient matrix
  assembly 
%
%  solve system  using Cholesky
f=a\f;
%
%  compute error
  err
end
 
