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
for i=1:1000
[xt, yt, zt]=target_motion_model2(dt,tempx,tempy,tempz);
tempx=xt;
tempy=yt;
tempz=zt;
X=[X xt];
Y=[Y yt];
Z=[Z zt];
S.x=xt+0.9*randn(1);
S.y=yt+0.9*randn(1);
S.z=zt+0.9*randn(1);
S_traj=[S_traj S];
plot(i,yt,'r.-');
hold on
plot(i,S.y,'b.-');
hold on
end
grid on
z_max=20;
z_min=-20;
x_max=20;
x_min=-20;
y_max=20;
y_min=-20;
%[St,Wt]=particle_filter_target(St_prior,Wt_prior,dt,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection);
% Xt(1)=S_traj(2).x;
% Xt(2)=S_traj(2).y;
% Xt(3)=S_traj(2).z;
Xt=[ 2.4431   -1.6725    0.7908];
% x_target=X(2);
% y_target=Y(2);
% z_target=Z(2);
 x_target= -2.1280;
 y_target=2.5223;
 z_target=-3.0416;
detect=true;
Pro_ztlxt=target_measurement_model(x_target,y_target,z_target,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detect)

% Xt(1)=S_traj(300).x;
% Xt(2)=S_traj(300).y;
% Xt(3)=S_traj(300).z;
% detect=true;
% Pro_ztlxt=target_measurement_model(x_target,y_target,z_target,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detect)
% 

