function [Pro_ztlxt]=needle_measurement_model(INT,x_int,y_int,x_max,x_min,y_max,y_min,Xt,detect)
%Xt=[x_image, y_image, z_image] is the position on the image of ultra sound
%imaging
alpha_d=0.9; % frequency percentage of detection
alpha_fp=0.1;%frequency percentage of false positive
%%
if INT==true
Pro_detect=alpha_d+alpha_fp;%Probability of zt1=true given Xt or Trustworthiness of the ultra sound detect the object
else
Pro_detect=alpha_fp;
end
Pro_nodetect=1-Pro_detect; %Probability of zt1=false given Xt. or the Trustworthiness of the ultra sound do not detect the object
%%
% target_x=mean(Xt(:,1));
% sigma_x=sqrt(var(Xt(:,1)));
sigma_x=1; %standard diviation of the noise in x_image
Pro_zt2_real_detect=normpdf(Xt(1),x_int,sigma_x);
Pro_zt2_false_detect=unifpdf(Xt(1),x_min,x_max);
%
if INT==true
Pro_real_detect=alph_d/(alpha_d+alpha_fp);
Pro_false_detect=alpha_fp/(alpha_d+alpha_fp);
else
Pro_real_detect=0;
Pro_false_detect=1;
end
% 
Pro_zt2=Pro_zt2_real_detect*Pro_real_detect+Pro_zt2_false_detect*Pro_false_detect;
%% 
sigma_y=1;%standard diviation of the noise in x_image
Pro_zt3_real_detect=normpdf(Xt(2),y_int,sigma_y);
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

