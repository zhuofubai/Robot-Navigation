function [xt, yt, zt]=target_motion_model2(dt,xt_prior, yt_prior, zt_prior)
%dt:  time step
%xt_prior, yt_prior, zt_prior : the target position at time t-1
%xt, yt, zt: the target position at time t
sigma=0.05; 
d_x=0.5;%velocity of target according to x coordinate
d_y=0.5;%velocity of target according to y coordinate
d_z=0.5;%velocity of target according to z coordinate
%%

z=randn(1);
index=randperm(2);
if index(1)==1
    delta_x=z*sigma+d_x*dt;
else 
    delta_x=-(z*sigma+d_x*dt);
end

z=randn(1);
index=randperm(2);
if index(1)==1
    delta_y=z*sigma+d_y*dt;
else 
    delta_y=-z*sigma-d_y*dt;
end

z=randn(1);
index=randperm(2);
if index(1)==1
    delta_z=z*sigma+d_z*dt;
else 
    delta_z=-z*sigma-d_z*dt;
end

xt=xt_prior+delta_x;
yt=yt_prior+delta_y;
zt=zt_prior+delta_z;

% xt=xt_prior+randn(1)*0.5;
% yt=yt_prior+randn(1)*0.5;
% zt=zt_prior+randn(1)*0.5;