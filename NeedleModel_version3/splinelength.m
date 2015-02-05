function len = splinelength(ABC,Xp1,Xp2)
a0 = ABC(1,1);
a1 = ABC(1,2);
a2 = ABC(1,3);
b0 = ABC(2,1);
b1 = ABC(2,2);
b2 = ABC(2,3);
c0 = ABC(3,1);
c1 = ABC(3,2);
c2 = ABC(3,3);

t11=Xp1(1);
t12=Xp1(2);
t13=Xp1(3);

t21=Xp2(1);
t22=Xp2(2);
t23=Xp2(3);

%find the length of the spline
%to avoid singularities, do it in rho-mu-tangent parameterization
% a1 -> rho*t11
% b1 -> rho*t12
% c1 -> rho*t13
% a2 -> (mu*t21 - rho*t11)/2
% b2 -> (mu*t22 - rho*t12)/2
% c2 -> (mu*t23 - rho*t13)/2

% rhos = [];
% if t11~=0
%     rhos =[rhos a1/t11];
% end
% if t12~=0
%     rhos = [rhos b1/t12];
% end
% if t13~=0
%     rhos = [rhos c1/t13];
% end
% if isempty(rhos)
%     rho = NaN;
% else
%     rho = max(rhos);
% end
% 
% mus = [];
% if t21 ~= 0
%     mus = [mus 2*a2+rho*t11/t21];
% end
% if t22 ~= 0
%     mus = [mus 2*b2+rho*t12/t22];
% end
% if t23 ~= 0
%     mus = [mus 2*c2+rho*t13/t23];
% end
% if isempty(mus)
%     mu = NaN;
% else
%     mu = max(mus);
% end



if t11~=0
    rho = a1/t11;
elseif t12~=0
    rho = b1/t12;
elseif t13~=0
    rho = c1/t13;
else
    rho = NaN;
end

if t21 ~= 0
    mu = 2*a2+rho*t11/t21;
elseif t22 ~= 0
    mu = 2*b2+rho*t12/t22;
elseif t23 ~= 0
    mu = 2*c2+rho*t13/t23;
else
    mu = NaN;
end

EPSILON = 1e-8;
if t11 < EPSILON
    t11=EPSILON;
end
if t12 < EPSILON
    t12=EPSILON;
end
if t13 < EPSILON
    t13=EPSILON;
end
if t21 < EPSILON
    t21=EPSILON;
end
if t22 < EPSILON
    t22=EPSILON;
end
if t23 < EPSILON
    t23=EPSILON;
end


len = ((rho*(rho - mu*(t11*t21 + t12*t22 + t13*t23)))/(mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23)) + ...
  (mu^2*(mu - rho*(t11*t21 + t12*t22 + t13*t23)))/(mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23)) + ...
  (mu^2*rho^2*(t13^2*(t21^2 + t22^2) - 2*t11*t13*t21*t23 - 2*t12*t22*(t11*t21 + t13*t23) + t12^2*(t21^2 + t23^2) + ...
     t11^2*(t22^2 + t23^2))*log(mu^2 - mu*rho*(t11*t21 + t12*t22 + t13*t23) + ...
      mu*sqrt(mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23))))/...
   (mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23))^(3/2) - ...
  (mu^2*rho^2*(t13^2*(t21^2 + t22^2) - 2*t11*t13*t21*t23 - 2*t12*t22*(t11*t21 + t13*t23) + t12^2*(t21^2 + t23^2) + ...
     t11^2*(t22^2 + t23^2))*log(-rho^2 + mu*rho*(t11*t21 + t12*t22 + t13*t23) + ...
      rho*sqrt(mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23))))/...
   (mu^2 + rho^2 - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23))^(3/2))/2;

if isnan(len)
    len = 0;
end
