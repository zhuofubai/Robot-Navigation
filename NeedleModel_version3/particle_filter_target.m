function [St,Wt]=particle_filter_needle(St_prior,Wt_prior,dt,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection)
St=St_prior;
Wt=Wt_prior;
% St_prior=low_variance_sampler(St_prior,Wt_prior);
len=length(Wt_prior);
weight_sum=0;
for i=1:len
  [St(i).x, St(i).y, St(i).z]=target_motion_model2(dt,St_prior(i).x, St_prior(i).y, St_prior(i).z);
 
  Wt(i)= target_measurement_model(St(i).x,St(i).y,St(i).z,z_max,z_min,x_max,x_min,y_max,y_min,Xt,detection(i));
  weight_sum=weight_sum+Wt(i);
end
Wt=Wt/weight_sum;

