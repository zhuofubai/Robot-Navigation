clear all
close all
%set up an initially straight needle pointing straight along the z-axis
clear all
% close all
clf

%set up important variables
doVisualization = 1;
nControlPoints = 5;
nSteps = 20;

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
%declare some globals for the nonlinear constraint function to use
for iStep=1:20
[G,Gp,ABC]=needle_motion_model_3(G, Gp, ABC,iStep,doVisualization);
end
z_target=G(3,3);
[xint yint INT]=find_position(z_target, ABC, G);