close all
clear all

bounds = [0 100 0 100]; % x and y min and max

Pneedle_init = [0 0]; % starting needle position
Rneedle_init = [0 1]; % starting needle direction

Pneedle = Pneedle_init;
Rneedle = Rneedle_init;

goalX = [100 100]; % target position



%centers and radii for the (circular) targets 
objX = [25 50;
        50 25];
objRad = [10 10]';
nObj = length(objRad);

nPts = 1000;
objBounds = zeros(nObj,2,nPts);

for i=1:nObj
    for j=1:nPts
        theta = 2*pi*(j-1) / nPts;
        objBounds(i,:,j) = objX(i,:) + (objRad .* [cos(theta);sin(theta)])';
    end
end

for i=1:nObj
    plot(objBounds(i,:,:));
    hold on
end
