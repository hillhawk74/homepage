% M-file for computing the L2 and H1 errors
% using 13 point quadrature rule
%
el2=0;
eh1=0;
for it=1:nel
   quad13
   for iq=1:13
      ar=area(it)*we(iq);
     	x=xe( iq);
     	y=ye( iq);
     	uh=0;
     	uhx=0;
     	uhy=0;
     	for in=1:nlocn
	        ip=node(it,in);
	        [bb,bx,by]=qbf(x,y,it,in,xc,yc,node);
	        i=indx(ip);
	        if i> 0
	            uh=uh+bb*f(i);
		        uhx=uhx+bx*f(i);
		        uhy=uhy+by*f(i);
           else
             uh=uh+bb*ubdry(xc(ip),yc(ip));
             uhx=uhx+bx*ubdry(xc(ip),yc(ip));
             uhy=uhy+by*ubdry(xc(ip),yc(ip));
	       end
    	end
    	[uex,uexx,uexy]=exact(x,y);
	    el2=el2+ar*(uh-uex)^2;
    	eh1=eh1+ar*( (uhx-uexx)^2+(uhy-uexy)^2 );
   end
end
el2=sqrt(el2);
eh1=sqrt(eh1);
disp(' ')
disp(sprintf('  l2-error is  %12.5e',el2))
disp(sprintf('  h1-error is  %12.5e',eh1))
