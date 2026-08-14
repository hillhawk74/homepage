% M-file to assemble matrix using the midpoint rule for quadrature
%
% 
a=zeros(nunk,nunk);
f=zeros(nunk,1); 
for it=1:nel
   for iq=1:nq
      x=xq(it,iq);
	     y=yq(it,iq);
	     ar=area(it)*wq(iq);
	     for in=1:nlocn
	        ip=node(it,in);
	        i=indx(ip);
	        if i>0
	          [bb,bx,by]=qbf(x,y,it,in,xc,yc,node);
		         [rhs]=rhsfun(x,y);
 		        f(i)=f(i)+rhs*bb*ar;
		         for inn=1:nlocn
		            ipp=node(it,inn);
			           [bbb,bbx,bby]=qbf(x,y,it,inn,xc,yc,node);
	     	      	j=indx(ipp); 
			           aij=bx*bbx+by*bby;
              if j>0
			              a(i,j)=a(i,j)+aij*ar; 
               else 
                 f(i)=f(i)-aij*ubdry(xc(ipp),yc(ipp))*ar ;
			          end
		         end
         end
      end
   end
end 
%
  

