function [St2]=low_variance_sampler(St,Wt)
    St2=[];
    M=length(Wt);
    r=(1/M)*rand(1);
    c=Wt(1);
    i=1;
    for m=1 : M
        U=r+(m-1)*(1/M);
        while U>c
            i=i+1;
            c=c+Wt(i);
        end
        St2=[St2 St(i)];
    end