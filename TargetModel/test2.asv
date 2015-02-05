close all
clear all
S=struct('x',[],'y',[],'z',[]);
S_traj=[];
X=[];
Y=[];
Z=[];
x=0;
y=0;
z=0;
tempx=0;
tempy=0;
tempz=0;
dt=1;
S_image=[];
ParticleSize=1000;
for i=1:10
    [xt, yt, zt]=target_motion_model2(dt,tempx,tempy,tempz);
    tempx=xt;
    tempy=yt;
    tempz=zt;
    X=[X xt];
    Y=[Y yt];
    Z=[Z zt];
    plot(i,yt,'r.-');
    hold on
    Zt.x=xt+0.01*randn(1);
    Zt.y=yt+0.01*randn(1);
    Zt.z=zt+0.01*randn(1);
    S_image=[S_image Zt];
end


z_max=5;
z_min=-5;
x_max=5;
x_min=-5;
y_max=5;
y_min=-5;
St_prior=[];
Wt_prior=zeros(1,ParticleSize);
detection=zeros(1,ParticleSize);

for i=1:ParticleSize
%     S.x=randn(1)*0.1;
%     S.y=randn(1)*0.1;
%     S.z=randn(1)*0.1;
    mms=randn(1,3)*0.1;
    S.x=mms(1); S.y=mms(2); S.z=mms(3);
    St_prior=[St_prior S];
    Wt_prior(i)=1/ParticleSize;
    detection(i)=true;
%     p=randperm(100);
%     if(p(1)==1)
%         detection(i)=false;
%     else
%         detection(i)=true;
%     end
end
%%
% Xt(1)=S_image(2).x;
% Xt(2)=S_image(2).y;
% Xt(3)=S_image(2).z;
% [St,Wt]=particle_filter_target(S_prior,Wt_prior,dt,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection);
% real=[X(2) Y(2) Z(2)]
% maxprob_index=find(Wt==max(Wt));
% maxprob=St(maxprob_index)
% St=low_variance_sampler(St,Wt);
% xsize=20;ysize=20;zsize=20;
% [density A]=discrete_histgram(xsize,ysize,zsize, St);
% density
% Xt
% %%
% Xt(1)=S_image(3).x;
% Xt(2)=S_image(3).y;
% Xt(3)=S_image(3).z;
% [St,Wt]=particle_filter_target(St,Wt,dt,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection);
% real2=[X(3) Y(3) Z(3)]
% maxprob_index=find(Wt==max(Wt));
% maxprob2=St(maxprob_index)
% St=low_variance_sampler(St,Wt);
% xsize=20;ysize=20;zsize=20;
% [density2 A]=discrete_histgram(xsize,ysize,zsize, St);
% density2
% Xt
%%
St=St_prior;
Wt=Wt_prior;
density=[];
real=[];
maxprob=[];
Parti_x=zeros(ParticleSize,1);
Parti_y=zeros(ParticleSize,1);
Parti_z=zeros(ParticleSize,1);
figure,
for i=1:length(S_image)
Xt(1)=S_image(i).x;
Xt(2)=S_image(i).y;
Xt(3)=S_image(i).z;
[St,Wt]=particle_filter_target(St,Wt,dt,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection);
real2=[X(i) Y(i) Z(i)];
real=[real;real2];
maxprob_index=find(Wt==max(Wt));
maxprob2=St(maxprob_index);
maxprob=[maxprob maxprob2];
St=low_variance_sampler(St,Wt);
%xsize=20;ysize=20;zsize=20;
%[density2 A]=discrete_histgram(xsize,ysize,zsize, St,Wt);
%density=[density density2];
% %
for j=1:ParticleSize
Parti_x(j)=St(j).x;
Parti_y(j)=St(j).y;
Parti_z(j)=St(j).z;
end
% density3(i,1)=density2.x; density3(i,2)=density2.y; density3(i,3)=density2.z;
[fx,xi] = ksdensity(Parti_x);
[fy,yi] = ksdensity(Parti_y);
[fz,zi] = ksdensity(Parti_z);
x_index=find(fx==max(fx));
y_index=find(fy==max(fy));
z_index=find(fz==max(fz));
density3(i,1)=xi(x_index); density3(i,2)=yi(y_index); density3(i,3)=zi(z_index);
 maxprob3(i,1)=maxprob2.x; maxprob3(i,2)=maxprob2.y; maxprob3(i,3)=maxprob2.z;
 real3(i,1)=real2(1); real3(i,2)=real2(2); real3(i,3)=real2(3);
 

figure,
for j=1:ParticleSize
     plot3(St(j).x,St(j).y,St(j).z,'bo');
     hold on
end
 grid on
 
% plot3 (density2.x, density2.y, density2.z,'ro');
% hold on
% plot3 (maxprob2.x, maxprob2.y, maxprob2.z,'go');
% hold on
% plot3 (real2(1), real2(2), real2(3),'bo');
% hold on
end
figure,
plot3(density3(:,1),density3(:,2),density3(:,3),'bo-');
hold on
plot3(maxprob3(:,1),maxprob3(:,2),maxprob3(:,3),'go-');
hold on
plot3(real3(:,1),real3(:,2),real3(:,3),'ro-');
grid on
%grid on
%%
% figure,
% for i=1:ParticleSize
%      plot3(St(i).x,St(i).y,St(i).z,'bo');
%      hold on
% end
%  grid on