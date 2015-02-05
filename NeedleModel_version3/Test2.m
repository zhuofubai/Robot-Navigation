%Test2 Particle Filter for Needle Motion model
clear all
close all
%set up an initially straight needle pointing straight along the z-axis
clear all
% close all
clf

%set up important variables
doVisualization = 1;
nControlPoints = 3;
nSteps = 5;

Ginit = zeros(nControlPoints,3);
Ginit(:,3) = 0:1:nControlPoints-1;

Gpinit = zeros(nControlPoints,3);
Gpinit(:,3) = 1;

%set the current needle position and tangents
G = Ginit;
Gp = Gpinit;

%set up the initial parameters for the splines
%intially, they are all equal, and are straight lines pointing in the z
%direction
ABC = zeros(3,3,nControlPoints-1);
for i=1:nControlPoints-1
    ABC(:,:,i) = [0 0 0;
        0 0 0;
        i-1 1 0];
end
ABCinit = ABC;

global totalDelta 
totalDelta=0;

%acatual needle shape in 10 timestep
%nSteps=10;
temp=struct('G',[],'Gp',[],'ABC',[]);
actual_needle_motion=[];
Z_img=[];
X_img=[];
Y_img=[];
xx=0;yy=0;zz=0;
for iStep=1:nSteps
    
    nSteps=5;
    nControlPoints=3;
    doVisualization=1;
[G,Gp,ABC]=needle_motion_model_3(nSteps,nControlPoints,G, Gp, ABC,iStep,doVisualization);
temp.G=G;
temp.Gp=Gp;
temp.ABC=ABC;
actual_needle_motion=[actual_needle_motion temp];
zz=(G(nControlPoints,3)+G(nControlPoints-1,3))/2;
Z_img=[Z_img zz];
[xx yy INT]=find_position(zz,ABC, G);
xx=xx+randn(1)*0.02;
yy=yy+randn(1)*0.02;
X_img=[X_img xx];
Y_img=[Y_img yy];
end
%%
%initialize particles
hold on
Needle_set=[];

doVisualization=0;
ParticleNum=100*nControlPoints;
Wt=zeros(1,ParticleNum);
for i=1:ParticleNum
    temp.G=Ginit;
    temp.Gp=Gpinit;
    temp.ABC=ABCinit;
    Needle_set=[Needle_set temp];
    Wt(i)=1/ParticleNum;
end
%%
%test the performance of particle filters
x_max=5; x_min=-5;y_max=5;y_min=-5;
iStep=1;
weight_sum=0;
X_set=zeros(1,ParticleNum);Y_set=zeros(1,ParticleNum);
Xt(1)=X_img(1);Xt(2)=Y_img(2);Xt(3)=Z_img(1);
z_target=Z_img(1);
figure,
for i=1:length(Needle_set)
    nSteps=10;
    nControlPoints=3;
    doVisualization=0;
    [G,Gp,ABC]=needle_motion_model_3(nSteps,nControlPoints,Needle_set(i).G, Needle_set(i).Gp, Needle_set(i).ABC,iStep,doVisualization);
    Needle_set(i).G=G;
    Needle_set(i).Gp=Gp; 
    Needle_set(i).ABC=ABC;   
    [xint yint INT]=find_position(z_target, Needle_set(i).ABC, Needle_set(i).G);
    detect=true;
    Wt(i)=needle_measurement_model(INT,xint,yint,x_max,x_min,y_max,y_min,Xt,detect);
    weight_sum=weight_sum+Wt(i);
    X_set(i)= xint; Y_set(i)=yint;
%     i   
end
%%
Wt=Wt/weight_sum;
[max_val max_ind]=max(Wt);
max_x = X_set(max_ind); max_y=Y_set(max_ind); 
Xt
[max_x max_y]
[Needle_set]=low_variance_sampler(Needle_set,Wt);

[fx,xi] = ksdensity(X_set);
[fy,yi] = ksdensity(Y_set);
[max_x_den, max_x_den_ind]=max(fx);
[max_y_den, max_y_den_ind]=max(fy);
xi(max_x_den_ind)
yi(max_y_den_ind)

%%
%plot 
        XData = [];
        YData = [];
        ZData = [];
        ABC_p=actual_needle_motion(1).ABC;
        G_p=actual_needle_motion(1).G;
        for i=1:nControlPoints-1
            a0 = ABC_p(1,1,i);
            a1 = ABC_p(1,2,i);
            a2 = ABC_p(1,3,i);
            b0 = ABC_p(2,1,i);
            b1 = ABC_p(2,2,i);
            b2 = ABC_p(2,3,i);
            c0 = ABC_p(3,1,i);
            c1 = ABC_p(3,2,i);
            c2 = ABC_p(3,3,i);
            l = 0:.01:1;
            XData = [XData a0 + a1*l + a2*l.^2];
            YData = [YData b0 + b1*l + b2*l.^2];
            ZData = [ZData c0 + c1*l + c2*l.^2];
        end
        plot3(G_p(:,3),G_p(:,1),G_p(:,2),'bo');
        plot3(ZData,XData,YData,'b');
        hold on
for i=1:length(Needle_set)
    plot3(Needle_set(i).G(:,3),Needle_set(i).G(:,1),Needle_set(i).G(:,2),'ro');
    hold on
    grid on
end