function [V_x V_y V_theta]=unicylce_control(xd, yd,theta)
kp=1;
ktheta=1;
%theta_d=atan2{xd,yd}-theta;

u1=kp*(xd*cos(theta)+y*sin(theta));
if(xd==0&&yd==0)
u2=theta;
else
u2=ktheta*atan2(yd,xd)-theta;
end
V_x=cos(theta)*u1;
V_y=sin(theta)*u1;
V_theta=u2;
