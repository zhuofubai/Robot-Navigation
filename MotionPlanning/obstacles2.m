% close all
clf
clear all
close all
%% parameters for attraction/repulsion
Qstar = 10;
dstar = 5;
zeta = 1;
eta = 1000;
EPSILON = 1e-3;

mass = 1;
vel = [0 0];
ts = .5;
endtime = 1000;

%% set up targets, obstacles, etc


Pneedle_init = [0 0]; % starting needle position
Rneedle_init = [0 1]; % starting needle direction

Pneedle = Pneedle_init;
Rneedle = Rneedle_init;

goalX = [100 100]; % target position

%centers and radii for the (circular) obstacles
obsX = [55 50;
    20 25
    45 25];
obsRad = [10 10 10]';
nObs = length(obsRad);

nPts = 1000;
obsBounds1 = zeros(nObs,2,nPts);
obsBounds2 = zeros(nObs,2,nPts);

%parameters of nonholonomic constraints
kp=1;
ktheta=0.2;
theta2=pi/16;
%running the needle
for i=1:nObs
    for j=1:nPts
        theta = 2*pi*(j-1) / nPts;
%         obsBounds1(i,:,j) = obsX(i,:) + (obsRad .* [cos(theta);sin(theta)])';
        obsBounds1(i,1,j) = obsX(i,1) + obsRad(i) * cos(theta);
        obsBounds1(i,2,j) = obsX(i,2) + obsRad(i) * sin(theta);
%         obsBounds2(i,:,j) = obsX(i,:) + ((obsRad+Qstar) .* [cos(theta);sin(theta)])';
        obsBounds2(i,1,j) = obsX(i,1) + (obsRad(i)+Qstar) * cos(theta);
        obsBounds2(i,2,j) = obsX(i,2) + (obsRad(i)+Qstar) * sin(theta);
    end
end

for i=1:nObs
    hold on
    obs1 = squeeze(obsBounds1(i,:,:));
    plot(obs1(1,:),obs1(2,:),'b')
    obs2 = squeeze(obsBounds2(i,:,:));
    plot(obs2(1,:),obs2(2,:),'--b')
end

plot(goalX(1),goalX(2),'bo');
plot(Pneedle(1),Pneedle(2),'r.');
% legend('needle position','target',2);
% axis([-5 105 -5 105]);

for itime=1:endtime
    
    %% attractive force due to the target
    
    % find the vector that points to the target
    v = goalX-Pneedle;
    % find the distance to the target
    d = norm(v);
    
    if d < EPSILON
        break;
    end
    
    % now normalize to get a unit vector
    % v = v/norm(v);
    
    % find the attractive force
    
    if d <= dstar
        Fa = zeta * v; % quadratic potential
    else
        Fa = (dstar * zeta * v) / d;
    end
    
    %% repulsive force due to the obstacles
    Fr = zeros(2,nObs);
    for i=1:nObs
        
        % find the vector that points to the obstacle
        v = Pneedle-obsX(i,:);
        % find the distance to the target
        d = norm(v);
        % now normalize to get a unit vector
        %    v = v/norm(v);
        
        % the force will be proportional to the distance and in the direction
        % of the vector
        
        if d > Qstar+obsRad(i)
            f = 0;
        elseif d < obsRad(i)
            f = 1;
        else
            f = eta * (1/Qstar - 1/d) * d^-2;
        end
        Fr(:,i) = (v * f)';
    end
    
    %calculate acceleration
    Ftotal = Fa';
    for i=1:nObs
        Ftotal = Ftotal+Fr(:,i);
    end
    Ftotal = Ftotal';
    
    a = Ftotal/mass;
    vel = .5*a*ts^2;
   % nonholonomic contraints  
    if (vel(2)==0&&vel(1)==0)
    theta_d=theta2;
    else
    theta_d=atan2(vel(2),vel(1))-theta2;
    end
    u1=kp*[cos(theta2) sin(theta2)]*vel';
    u2=ktheta*theta_d;
    vel2=[cos(theta2)*u1 sin(theta2)*u1];
    
    
    Pneedle = Pneedle + vel2*ts;
    theta2=theta2+u2*ts;
    
    plot(Pneedle(1),Pneedle(2),'r.');
    title(sprintf('timestep %d, F = [%.2f %.2f]',itime, Ftotal(1),Ftotal(2)));
    pause(.001);
    
end