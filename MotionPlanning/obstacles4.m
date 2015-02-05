% close all
clf
clear all
close all


hold on
grid on

%% parameters for attraction/repulsion
Qstar = 20;
dstar = 5;
zeta = 1;
eta = 2500;
EPSILON = 1e0;
%dd
mass = 1;
vel = [0 0 0];
ts = .5;
endtime = 1000;

%% set up targets, obstacles, etc


Pneedle_init = [0 0 0]; % starting needle position
Rneedle_init = [0 0 1]; % starting needle direction
    
Pneedle = Pneedle_init;
Rneedle = Rneedle_init;
Rneedle = Rneedle/norm(Rneedle);

goalX = [100 100 100]; % target position

%centers and radii for the (circular) obstacles
obsX = [65 35 35;
    20 35 35;
    80 70 70];
obsRad = [10 10 10]';
nObs = length(obsRad);

% nPts = 1000;
% obsBounds1 = zeros(nObs,2,nPts);
% obsBounds2 = zeros(nObs,2,nPts);

%parameters of nonholonomic constraints
kp=2;
ktheta=0.05; %this controls the turning speed
kphi=0.05;
% generate obstacles
% for i=1:nObs
%     for j=1:nPts
%         theta = 2*pi*(j-1) / nPts;
%         %         obsBounds1(i,:,j) = obsX(i,:) + (obsRad .* [cos(theta);sin(theta)])';
%         obsBounds1(i,1,j) = obsX(i,1) + obsRad(i) * cos(theta);
%         obsBounds1(i,2,j) = obsX(i,2) + obsRad(i) * sin(theta);
%         %         obsBounds2(i,:,j) = obsX(i,:) + ((obsRad+Qstar) .* [cos(theta);sin(theta)])';
%         obsBounds2(i,1,j) = obsX(i,1) + (obsRad(i)+Qstar) * cos(theta);
%         obsBounds2(i,2,j) = obsX(i,2) + (obsRad(i)+Qstar) * sin(theta);
%     end
% end

%plot obstacles
% for i=1:nObs
%     hold on
%     obs1 = squeeze(obsBounds1(i,:,:));
%     plot(obs1(1,:),obs1(2,:),'b')
%     obs2 = squeeze(obsBounds2(i,:,:));
%     plot(obs2(1,:),obs2(2,:),'--b')
% end



for i=1:nObs
    [x,y,z] = sphere;
    x = x*obsRad(i);
    y = y*obsRad(i);
    z = z*obsRad(i);
    
    x = x+obsX(i,1);
    y = y+obsX(i,2);
    z = z+obsX(i,3);
    obsBounds1(:,:,i) = [x,y,z];
end

for i=1:nObs
    d = size(obsBounds1(:,:,i),1);
    
    x = obsBounds1(:,1:d,i);
    y = obsBounds1(:,d+1:2*d,i);
    z = obsBounds1(:,2*d+1:end,i);
    surf(x,y,z);
end

    
plot3(goalX(1),goalX(2),goalX(3),'bo');
plot3(Pneedle(1),Pneedle(2),Pneedle(3),'r.');
daspect([1 1 1]);
% legend('needle position','target',2);
axis([-5 105 -5 105 -5 105]);

    
HIT = false;
    
for itime=1:endtime
    Pneedlelast = Pneedle;
    
    %% attractive force due to the target
    
    % find the vector that points to the target
    v = goalX-Pneedle;
    % find the distance to the target
    d = norm(v);
    
    if d < EPSILON
        HIT = true;
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
    Fr = zeros(3,nObs);
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
        Ftotal = Ftotal-Fr(:,i);
    end
    Ftotal = Ftotal';
    
    a = Ftotal/mass;
    vel = .5*a*ts^2;
    %% nonholonomic contraints
    
    theta2 = real(atan2(Rneedle(2),Rneedle(1)) - atan2(0, 1)); %e.g. acos( needle vector dot xaxis )
    phi=real(atan2(Rneedle(3),sqrt(Rneedle(1)^2+Rneedle(2)^2))-atan2(0,1));
    %calculate u1
    u1=kp*[cos(theta2)*cos(phi) sin(theta2)*cos(phi) sin(phi)]*vel';
    %calculate u2 
    if (vel(2)==0&&vel(1)==0)
         theta_d=theta2;
     else
         theta_d=atan2(vel(2),vel(1))-theta2;
     end
    u2=ktheta*theta_d;
    %calculate u3
    u3=u1*cos(phi);
    %calculate u4
    vel_plannar=sqrt(vel(1)^2+vel(2)^2);
    if (vel(3)==0&&vel_plannar==0)
         phi_d=phi;
     else
         phi_d=atan2(vel(3),vel_plannar)-phi;
    end
     u4=kphi*phi_d;
    
     vel2=[cos(theta2)*u1*cos(phi) sin(theta2)*u1*cos(phi) u1*sin(phi)];
    
%vel2 = vel;
    
    Pneedle = Pneedle + vel2*ts;
    
    %if we're stopped, exit
%     if norm(Pneedle-Pneedlelast) < .0001
%         break;
%     end
    
     theta2=theta2+u2*ts;
     phi=phi+u4*ts;
     Rneedle = [cos(theta2)*cos(phi) sin(theta2)*cos(phi) sin(phi)];
    Rneedle = Rneedle/norm(Rneedle);
%%    
    plot3(Pneedle(1),Pneedle(2),Pneedle(3),'r.');
    title(sprintf('timestep %d, F = [%.2f %.2f]',itime, Ftotal(1),Ftotal(2)));
    pause(.005);
    
end
 plot3(Pneedle(1),Pneedle(2),Pneedle(3),'r.');
 if HIT
     title(sprintf('hit target at timestep %d, F = [%.2f %.2f]',itime, Ftotal(1),Ftotal(2)));
 else
     title(sprintf('stopped at timestep %d, F = [%.2f %.2f]',itime, Ftotal(1),Ftotal(2)));
 end