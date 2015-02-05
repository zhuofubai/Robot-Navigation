function [density A]=discrete_histgram(xsize,ysize,zsize, St)
bin_total=xsize*ysize*zsize;
A=zeros(xsize,ysize,zsize);
x=0;y=0;z=0;
for i=1:length(St)
    if (St(i).x>=0)
      x=fix(St(i).x)+1+xsize/2;
    else
       x=fix(St(i).x)+xsize/2; 
    end
    
    if (St(i).y>=0)
      y=fix(St(i).y)+1+ysize/2;
    else
       y=fix(St(i).y)+ysize/2; 
    end
    
    if (St(i).z>=0)
      z=fix(St(i).z)+1+zsize/2;
    else
      z=fix(St(i).z)+zsize/2; 
    end
      A(x,y,z)=A(x,y,z)+1;
end
siz=[xsize ysize zsize];
ind=find(A==max(max(max(A))));
ind=ind(1);
[I,J,K] = ind2sub(siz,ind);
density.x=I-xsize/2;
density.y=J-ysize/2;
density.z=K-zsize/2;