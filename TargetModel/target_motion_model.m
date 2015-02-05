function [Pro_xt,Pro_yt,Pro_zt]=target_motion_model(dt,xt_prior, yt_prior, zt_prior, xt, yt, zt)
%dt:  time step
%xt_prior, yt_prior, zt_prior : the target position at time t-1
%xt, yt, zt: the target position at time t
sigma=1; 
d_x=1;%velocity of target according to x coordinate
d_y=1;%velocity of target according to y coordinate
d_z=1;%velocity of target according to z coordinate
%%

delta_x=d_x*dt+xt_prior;
K=sigma*dt;
Pro_xt=exp((-1)*(xt-delta_x)^2/(2*K))/sqrt(2*pi*K);

delta_y=d_y*dt+yt_prior;
K=sigma*dt;
Pro_yt=exp((-1)*(yt-delta_y)^2/(2*K))/sqrt(2*pi*K);

delta_z=d_z*dt+zt_prior;
K=sigma*dt;
Pro_zt=exp((-1)*(zt-delta_z)^2/(2*K))/sqrt(2*pi*K);