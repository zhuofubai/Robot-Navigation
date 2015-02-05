function [Pro_ztlxt]=target_measurement_model(x_target,y_target,z_target,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detect)
%Xt=[x_image, y_image, z_image] is the position on the image of ultra sound
%imaging
z_image=Xt(3);
epsilon1=0;  %¦Å1
epsilon2=1/(z_max-z_min);% ¦Å2
Area=intesection_area(z_target,z_image);%Area of the target in the image.
Pro_det=Detection_probability(Area,epsilon1);  % probability of detection
Pro_false_pos=1/(z_max-z_min); %  parobability of false positive
alpha_d=0.9; % frequency percentage of detection
alpha_fp=0.01;%frequency percentage of false positive
%%
Pro_detect=Pro_det*alpha_d+Pro_false_pos*alpha_fp;%Probability of zt1=true given Xt or Trustworthiness of the ultra sound detect the object
Pro_nodetect=1-Pro_detect; %Probability of zt1=false given Xt. or the Trustworthiness of the ultra sound do not detect the object
%%
% target_x=mean(Xt(:,1));
% sigma_x=sqrt(var(Xt(:,1)));
sigma_x=0.03; %standard diviation of the noise in x_image
Pro_zt2_real_detect=normpdf(Xt(1),x_target,sigma_x);
Pro_zt2_false_detect=unifpdf(Xt(1),x_min,x_max);
%
Pro_real_detect=alpha_d*Pro_det/(Pro_det*alpha_d+Pro_false_pos*alpha_fp);
Pro_false_detect=Pro_false_pos*alpha_fp/(Pro_det*alpha_d+Pro_false_pos*alpha_fp);
% 
Pro_zt2=Pro_zt2_real_detect*Pro_real_detect+Pro_zt2_false_detect*Pro_false_detect;
%% 
sigma_y=0.03;%standard diviation of the noise in x_image
Pro_zt3_real_detect=normpdf(Xt(2),y_target,sigma_y);
Pro_zt3_false_detect=unifpdf(Xt(2),y_min,y_max);
% 
Pro_zt3=Pro_zt3_real_detect*Pro_real_detect+Pro_zt3_false_detect*Pro_false_detect;
%% 

if detect==true
Pro_ztlxt=Pro_detect*Pro_zt2*Pro_zt3;
else
  Pro_ztlxt = Pro_nodetect;
end
%%
function [Area]=intesection_area(z_target,z_image)
R=3;%actual radius of the target
delta=abs(z_target-z_image);
if delta<R
    Area=pi*(R^2-abs(z_target-z_image)^2);
else
    Area=0;
end
function [P_det]=Detection_probability(Area,epsilon1)
R=3;
A0=pi*R*R;%actual area of the target
if Area<A0
    P_det=(1-epsilon1)*Area/A0;
else
    P_det=1;
end
