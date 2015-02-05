function [xd,yd]=potential_fun(Pneedle, goalX,obsX,obsRad)
%% attractive force due to the target
Qstar = 5;
dstar = 5;
zeta = 1;
eta = 10000;
% find the vector that points to the target
v = Pneedle-goalX;
% find the distance to the target
d = norm(v);
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
%     v = v/norm(v);
    
    % the force will be proportional to the distance and in the direction
    % of the vector
 
    if d > Qstar+obsRad(i)
        f = 0;
    elseif d < obsRad(i)
        f = 100; 
    else        
        f = eta * (1/Qstar - 1/d) * d^-2;
    end
    Fr(i,:) = v * f;    %v=v/norm(v)
end
%% vortex field
