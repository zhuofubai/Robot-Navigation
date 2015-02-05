%give the value of z cooridate, compute the corresponding x and y coordinate
function [xint yint INT]=find_position(z_target, ABC, G)
INT=true;
xint=0;yint=0;
index=find(G(:,3)>z_target);
TF = isempty(index);
if(TF==1)
    if G(5,3)==z_target
    spline_index=4;
    else
    INT=false;
    end
else
spline_index=index(1)-1;
end
if (length(index)>1&&index(1)==index(2))    
        spline_index=index(2)-1;
end


if (spline_index==0)
    INT=false;
end
if INT==true
    cp=spline_index;
    ABCcur = ABC(:,:,cp);
    %advance the control point along the spline that starts at cp
    a0 = ABCcur(1,1);
    a1 = ABCcur(1,2);
    a2 = ABCcur(1,3);
    b0 = ABCcur(2,1);
    b1 = ABCcur(2,2);
    b2 = ABCcur(2,3);
    c0 = ABCcur(3,1);
    c1 = ABCcur(3,2);
    c2 = ABCcur(3,3);
    %         Gcur(1) = a0 + a1*delta + a2*delta^2;
    %         Gcur(2) = b0 + b1*delta + b2*delta^2;
    %         Gcur(3) = c0 + c1*delta + c2*delta^2;
    %delta=0;
    if c2==0
        if c1==0
            delta=0;
            
        else
            delta=(z_target-c0)/c1;
        end
    else
        p=[c2 c1 c0-z_target];
        ans=roots(p);
        ans2=ans(find((ans<1)&(ans>-1)));
        
        if length(ans2)==2
            delta=ans2(1);
        elseif length(ans2)==0
            INT=false;
            delta=0;
            
        else
            delta=ans2;
        end
    end
    if delta>1
    delta
    end
    xint = a0 + a1*delta + a2*delta^2;
    yint = b0 + b1*delta + b2*delta^2;
end

