c *****************************************************************
c
c  This program solves -(laplacian u) =f in a box in R2 where 
c  homogeneous Dirichlet boundary conditions are imposed.  The
c  code uses continuous piecewise quadratic basis functions on
c  triangles.
c
c *****************************************************************
c *****************************************************************
c
c  INPUT
c    nx    -  number of points in the x-direction
c    ny    -  number of points in the y-direction
c    xl,xr -  left and right x-coordinates of domain
c    yb,yt -  bottom and top y-coordinates of domain
c    nnodes-  number of nodes per element (for quadratics it's 6)
c    nq    -  number of quadrature points per element used in the assembly
c    nqe   -  number of quadrature points per element used in the error routine
c    lda   -  leading dimension of a(i,j) 
c    lde   -  leading dimension of arrays dimensioned by number of elements
c
c *****************************************************************
c
c  GEOMETRY
c    xc(i)  -  x-coordinate of global node number i
c    yc(i)  -  y-coordinate of global node number i
c    node(it,j)  - array giving global node number corresponding to
c                  local node number j of element number it
c    area(it)    - area of element it
c    indx(i)     - gives unknwon number at global node number i
c                  (i=0 if no unknown at node)
c    xq(it,i)    - x-coordinate of quadrature point j for element it
c    yq(it,i)    - y-coordinate of quadrature point j for element it
c
c    nunk        - number of unknowns
c    nel         - number of elements
c    np          - number of global nodes (or points)
c
c *****************************************************************
c
c  ASSEMBLY
c    a(i,j)  - coefficient matrix stored here using banded storage
c              dimensioned by the number of unknowns by total bandwidth
c    f(i)    - right hand side array which is overwritten by solver with
c              solution
c
c *****************************************************************
c *****************************************************************
c
      implicit double precision (a-h,o-z)
      dimension xc(1089),yc(1089),indx(1089),f(1089),wk(1089),
     .          a(961,129),
     .          area(512),xq(512,3),yq(512,3),node(512,6),wq(3),
     .          xe(512,13),ye(512,13),we(13)
c
c  input data
      lde=512
      lda=961
      nx=3
      ny=3
      do 500 ngrid=1,3
      if(ngrid.eq.1) nx=5
      if(ngrid.eq.2) nx=9 
	  if(ngrid.eq.3) nx=17
      ny=nx
      xl=0.
      xr=1.
      yb=0.
      yt=1.
      nq=3
      nqe=13
      nnodes=6
      write(0,1001)
      write(0,1002) nx,ny
c
c  set up geometry
      call geom(xc,yc,xq,yq,area,wq,node,indx,nel,nunk,np,nq,nnodes,
     .          xl,xr,yb,yt,nx,ny,lde)
c
c  assemble coefficient matrix and right-hand side
      call assem(a,f,ib,xc,yc,xq,yq,wq,area,node,indx,
     .           nnodes,nunk,nq,nel,lde,lda)
c
c  solve system using banded solver
      call bepp(lda,nunk,ib,ib,a,1,f,wk,ierr)
      if(ierr.eq.0) go to 100
      write(0,1005)
	  write(*,1005)
      stop
  100 continue
c
c  calculate error using 13 point quadrature rule
      call quad13(xc,yc,node,nel,lde,we,xe,ye)
      call eror(el2,eh1,we,xe,ye,area,node,indx,
     .          xc,yc,f,nel,nnodes,lde,nqe)
  500 continue
c
      stop
 1005 format('matrix is singular')
 1001 format(/,'solution of -uxx+uyy=f on unit box',/,
     .        '                   u=0 on boundary',/
     .        'piecewise quadratics on triangles are used',//)
 1002 format(/'number of points in x-direction =',i5,
     .       /'number of points in y-direction =',i5)
      end
c
c
c ****************************************************************
c *****************************************************************
      subroutine geom(xc,yc,xq,yq,area,wq,node,indx,nel,nunk,np,nq,
     .                nnodes,xl,xr,yb,yt,nx,ny,lde)
c  ***************************************************************
c   Sets up geometric information for a rectangular domain
c  ***************************************************************
c
      implicit double precision (a-h,o-z)
      dimension xc(1),yc(1),xq(lde,1),yq(lde,1),area(1),
     .          node(lde,1),indx(1),wq(1)
c
c  mapping functions for variable grid
      grdx(x)=x
      grdy(y)=y
c
c  set up grid information and unknown counter array
      nxm=nx-1
      hx=(xr-xl)/float(nxm)
      nym=ny-1
      hy=(yt-yb)/float(nym)
      yy=-hy/2.
      nx2=nx+nxm
      ny2=ny+nym
      iu=0
      ip=0
      it1=-1
      do 150 jc=1,ny2
        yy=yy+hy/2.
        xx=-hx/2.
        iy=mod(jc,2)
        do 140 ic=1,nx2
          xx=xx+hx/2.
          ix=mod(ic,2)
          ip=ip+1
          xc(ip)=grdx(xx)
          yc(ip)=grdy(yy)
          if(ix.eq.1.and.iy.eq.1) go to 110
          go to 120
  110     continue
          if(ic.eq.nx2.or.jc.eq.ny2) go to 120
          it1=it1+2
          it2=it1+1
          node(it1,1)=ip
          node(it1,2)=ip+nx2+nx2
          node(it1,3)=ip+nx2+nx2+2
          node(it1,4)=ip+nx2
          node(it1,5)=ip+nx2+nx2+1
          node(it1,6)=ip+nx2+1
          node(it2,1)=ip
          node(it2,2)=ip+2
          node(it2,3)=ip+nx2+nx2+2
          node(it2,4)=ip+1
          node(it2,5)=ip+nx2+2
          node(it2,6)=ip+nx2+1
 120      continue
            if(ic.eq.1.or.ic.eq.nx2) go to 130
            if(jc.eq.1.or.jc.eq.ny2) go to 130
          iu=iu+1
          indx(ip)=iu
          go to 140
 130      continue
          indx(ip)=0
 140    continue
 150  continue
      nunk=iu
      np=ip
      nel=it2
      write(0,1001) nunk,np,nel
c
c  set up quadrature information needed for assembly
       wq(1)=1./3.
       wq(2)=wq(1)
       wq(3)=wq(1)
       do 180 it=1,nel
         ip1=node(it,1)
         ip2=node(it,2)
         ip3=node(it,3)
         x1=xc(ip1)
         x2=xc(ip2)
         x3=xc(ip3)
         y1=yc(ip1)
         y2=yc(ip2)
         y3=yc(ip3)
         xq(it,1)=.5*(x1+x2)
         xq(it,2)=.5*(x2+x3)
         xq(it,3)=.5*(x1+x3)
         yq(it,1)=.5*(y1+y2)
         yq(it,2)=.5*(y2+y3)
         yq(it,3)=.5*(y1+y3)
         area(it)=.5*abs(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2))
  180  continue
c
 1001 format(/'number of unknowns=',i5,/'number of nodes=',i5,/,
     .       'number of elements=',i5)
       return
       end
c
c
c  ***************************************************************
c  ***************************************************************
      subroutine assem(a,f,ib,xc,yc,xq,yq,wq,area,node,indx,
     .           nnodes,nunk,nq,nel,lde,lda)
c  ***************************************************************
c
c  Sets up the coefficient matrix and right-hand side vector using
c  piecewise quadratics
c
c  ***************************************************************
c
      implicit double precision (a-h,o-z)
      dimension a(lda,1),f(1),xc(1),yc(1),xq(lde,1),yq(lde,1),
     .          area(1),wq(1),indx(1),node(lde,1)
c
c  determine bandwidth of coefficient matrix
      ib=0
      do 30 it=1,nel
        do 20 in=1,nnodes
          ip=node(it,in)
          i=indx(ip)
          if(i.eq.0) go to 20
          do 10 inn=1,nnodes
            ipp=node(it,inn)
            j=indx(ipp)
            ij=j-i
            if(ij.gt.ib) ib=ij
 10       continue
 20     continue
 30   continue
      ib1=ib+1
      ibt=ib+ib1
      write(0,1001) ibt
c
c  zero arrays
      do 40 i=1,nunk
        f(i)=0.
        do 40 j=1,ibt
          a(i,j)=0.0
 40   continue
c
c  assemble matrix and right hand side
      do 150 it=1,nel
        do 140 iq=1,nq
          x=xq(it,iq)
          y=yq(it,iq)
          ar=area(it)*wq(iq)
          do 130 in=1,nnodes
            ip=node(it,in)
            i=indx(ip)
            if(i.eq.0) go to 130
            call qbf(x,y,it,in,bb,bx,by,xc,yc,node,lde)
            ii=ib1-i
            f(i)=f(i)+rhs(x,y)*bb*ar
            do 120 inn=1,nnodes
              ipp=node(it,inn)
              call qbf(x,y,it,inn,bbb,bbx,bby,xc,yc,node,lde)
              j=indx(ipp)
               aij=bx*bbx+by*bby  
              if(j.eq.0) go to 115
              ij=ii+j
              a(i,ij)=a(i,ij)+aij*ar
              go to 120
 115          continue
              call exact(xc(ipp),yc(ipp),uex,uexx,uexy)
              f(i)=f(i)-aij*ar*uex
 120        continue
 130      continue
 140    continue
 150  continue
c
 1001 format(/'total bandwidth',i5)
      return
      end
c
c
c  ***************************************************************
c  ***************************************************************
      subroutine qbf(x,y,it,in,bb,bx,by,xc,yc,node,lde)
c  ***************************************************************
c  Quadratic basis functions
c  ***************************************************************
c
      implicit double precision (a-h,o-z)
      dimension node(lde,1),xc(1),yc(1)
      if(in.gt.3) go to 100
      in1=in
      in2=mod(in,3)+1
      in3=mod(in+1,3)+1
      i1=node(it,in1)
      i2=node(it,in2)
      i3=node(it,in3) 
      d=(xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
      t=1.+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d
      bb=t*(2.*t-1.)
      bx=(yc(i2)-yc(i3))*(4.*t-1.)/d
      by=(xc(i3)-xc(i2))*(4.*t-1.)/d
      return
 100  continue
      inn=in-3
      in1=inn
      in2=mod(inn,3)+1
      in3=mod(inn+1,3)+1
      i1=node(it,in1)
      i2=node(it,in2)
      i3=node(it,in3)
      j1=i2
      j2=i3
      j3=i1 
      d=(xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
      c=(xc(j2)-xc(j1))*(yc(j3)-yc(j1))-(xc(j3)-xc(j1))*(yc(j2)-yc(j1))  
      t=1.+((yc(i2)-yc(i3))*(x-xc(i1))+(xc(i3)-xc(i2))*(y-yc(i1)))/d
      s=1.+((yc(j2)-yc(j3))*(x-xc(j1))+(xc(j3)-xc(j2))*(y-yc(j1)))/c
      bb=4.*s*t
      bx=4.*(t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d)
      by=4.*(t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d)
      return
      end
c
c
c  ***************************************************************
c  ***************************************************************
      subroutine eror(el2,eh1,we,xe,ye,area,node,indx,
     .                xc,yc,f,nel,nnodes,lde,nq)
c  ***************************************************************
c
c  Calculates the L2 and H1-seminorm errors
c
c  ***************************************************************
c
      implicit double precision (a-h,o-z)
      dimension we(1),xe(lde,1),ye(lde,1),area(1),indx(1),f(1),
     .          xc(1),yc(1),node(lde,1)
      el2=0.
      eh1=0.
      do 150 it=1,nel
        do 140 iq=1,nq
          ar=area(it)*we(iq)
          x=xe(it,iq)
          y=ye(it,iq)
          uh=0.
          uhx=0.
          uhy=0.
          do 130 in=1,nnodes
            ip=node(it,in)
            call qbf(x,y,it,in,bb,bx,by,xc,yc,node,lde)
            x1=xc(ip)
            y1=yc(ip)
            i=indx(ip)
            if(i.le.0) go to 125
            uh=uh+bb*f(i)
            uhx=uhx+bx*f(i)
            uhy=uhy+by*f(i)
            go to 130
 125        continue
            call exact(x1,y1,uex,uexx,uexy)
            uh=uh+bb*uex
            uhx=uhx+bx*uex
            uhy=uhy+by*uex
 130      continue
          call exact(x,y,uex,uexx,uexy)
          el2=el2+(uh-uex)**2*ar
          eh1=eh1+((uhx-uexx)**2+(uhy-uexy)**2)*ar
 140    continue
 150  continue
      el2=sqrt(el2)
      eh1=sqrt(eh1)
      write(0,1001)
      write(0,1002) el2
      write(0,1003) eh1
	  write(*,1002) el2
      write(*,1003) eh1
      write(0,1001)
 1001 format(/,'*********************************************')
 1002 format(/,'L2-ERROR is  ',e12.5)
 1003 format('H1-SEMINORM ERROR is  ',e12.5)
      return 
      end
c
c
c  ***************************************************************
c  ***************************************************************
      subroutine quad13(xc,yc,node,nel,lde,we,xe,ye)
c  ***************************************************************
c
c  Sets up quadrature information for 13-point rule
c
c  ***************************************************************
      implicit double precision (a-h,o-z)
      dimension xc(1),yc(1),node(lde,1),we(13),xe(lde,1),ye(lde,1)
      do 10 i=1,3
        we(i)=.175615257433204
        ii=i+3
        we(ii)=.053347235608839
        ii=i+6
        iii=ii+3
        we(ii)=.077113760890257
        we(iii)=we(ii)
 10     continue
        we(13)=-.14957004446767
        z1=.479308067841923
        z2=.260345966079038
        z3=.869739794195568
        z4=.065130102902216
        z5=.638444188569809
        z6=.312865496004875
        z7=.048690315425316
        do 100 it=1,nel
          ip1=node(it,1)
          ip2=node(it,2)
          ip3=node(it,3)
          x1=xc(ip1)
          x2=xc(ip2)
          x3=xc(ip3)
          y1=yc(ip1)
          y2=yc(ip2)
          y3=yc(ip3)
          xe(it,1)=z1*x1+z2*x2+z2*x3
          ye(it,1)=z1*y1+z2*y2+z2*y3
          xe(it,2)=z2*x1+z1*x2+z2*x3
          ye(it,2)=z2*y1+z1*y2+z2*y3
          xe(it,3)=z2*x1+z2*x2+z1*x3
          ye(it,3)=z2*y1+z2*y2+z1*y3
          xe(it,4)=z3*x1+z4*x2+z4*x3
          ye(it,4)=z3*y1+z4*y2+z4*y3
          xe(it,5)=z4*x1+z3*x2+z4*x3
          ye(it,5)=z4*y1+z3*y2+z4*y3
          xe(it,6)=z4*x1+z4*x2+z3*x3
          ye(it,6)=z4*y1+z4*y2+z3*y3
          xe(it,7)=z5*x1+z6*x2+z7*x3
          ye(it,7)=z5*y1+z6*y2+z7*y3
          xe(it,8)=z5*x1+z7*x2+z6*x3
          ye(it,8)=z5*y1+z7*y2+z6*y3
          xe(it,9)=z6*x1+z5*x2+z7*x3
          ye(it,9)=z6*y1+z5*y2+z7*y3 
          xe(it,10)=z6*x1+z7*x2+z5*x3
          ye(it,10)=z6*y1+z7*y2+z5*y3
          xe(it,11)=z7*x1+z5*x2+z6*x3
          ye(it,11)=z7*y1+z5*y2+z6*y3
          xe(it,12)=z7*x1+z6*x2+z5*x3
          ye(it,12)=z7*y1+z6*y2+z5*y3 
          xe(it,13)=(x1+x2+x3)/3.
          ye(it,13)=(y1+y2+y3)/3.
 100      continue
          return
          end
c
c
c  ***************************************************************
c  ***************************************************************
      subroutine bepp(lda,n,m1,m2,a,nrhs,b,wk,ierr)
c  ***************************************************************
c
c  banded Gaussian elimination solver with partial pivoting
c
c  ***************************************************************
      implicit double precision (a-h,o-z)
      dimension a(lda,1),b(lda,1),wk(1)
      sicri=1.e-5
      ierr=0
      m=m1+m2+1
      do 60 i=1,n
        p=abs(a(i,1))
        do 50 j=2,m
          q=abs(a(i,j))
          if(q.gt.p)p=q
 50     continue
        fff=abs(p)
        if(fff.lt.sicri) go to 999
        wk(i)=p
 60   continue  
      do 80 i=1,m1
        je=m2+i
        js=m1+1-i
        do 70 j=1,je
          a(i,j)=a(i,js+j)
 70     continue
        je=je+1
        do 80 j=je,m
          a(i,j)=0.
 80     continue
        nm1=n-1
        l=m1
        do 140 k=1,nm1
          kp1=k+1
          if(l.lt.n)l=l+1
          i=k
          p=abs(a(k,1))/wk(k)
          do  90 j=kp1,l
            q=abs(a(j,1))/wk(j)
            if(q.le.p) go to 90
            p=q
            i=j
 90       continue
          fff=abs(p)
          if(fff.lt.sicri) go to 999
          if(i.eq.k) go to 120
          do 100 j=1,m
            x=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=x
 100      continue
          do 110 j=1,nrhs
            x=b(k,j)
            b(k,j)=b(i,j)
            b(i,j)=x
 110      continue
 120      continue
          do 135 i=kp1,l
            x=a(i,1)/a(k,1)
            do 130 j=2,m
              a(i,j-1)=a(i,j)-x*a(k,j)
 130        continue
            a(i,m)=0.
            do 135 j=1,nrhs
              b(i,j)=b(i,j)-x*b(k,j)
 135      continue
 140    continue
          fff=abs(a(n,1))
            if(fff.lt.sicri) go to 999
            np1=n+1
            do 180 j=1,nrhs
              l=1
              do 180 ii=1,n
                i=np1-ii
                fff=abs(a(i,1))
                if(fff.lt.sicri) go to 999
                im1=i-1
                x=b(i,j)
                if(i.eq.n) go to 170
                do 160 k=2,l
                  x=x-a(i,k)*b(im1+k,j)
 160            continue
 170            continue
                b(i,j)=x/a(i,1)
                if(l.lt.m) l=l+1
 180          continue
              go to 1000
 999          ierr=1
            write(*,1212) fff
 1212       format(' pivot element',e12.5)
 1000       continue
            return
            end
c  ***************************************************************
c  ***************************************************************
      subroutine exact(x,y,uex,uexx,uexy)
c  ***************************************************************
c
c  Calculates exact solution and its first derivatives
c  Must be supplied by user
c
c  ***************************************************************  
c
      implicit double precision (a-h,o-z)
      pi=4.*atan(1.)  
	  uex=(sin(pi*x)) *(sin(pi*y)) 
	  uexx= pi*cos(pi*x) *sin(pi*y)
	  uexy= pi*cos(pi*y)*sin(pi*x)  
      return
      end
c
c
c  ***************************************************************
c  ***************************************************************
      function rhs(x,y)
c  ***************************************************************
c
c  Gives the right-hand side of DE
c  Must be supplied by user
c
c  ***************************************************************
      implicit double precision (a-h,o-z) 
      pi=4.*atan(1.) 
      rhs=2.*pi*pi*sin(pi*x)*sin(pi*y) 
      return 
      end
c  ***************************************************************
 