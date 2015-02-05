function [c,ceq] = nonlconst(X)

%constraints:
%starting point: X1
%ending point: Xp2
%starting tangent: Xp1
%ending tangent: Xp2
%spline length: splinelen

global lambda1 lambda2 t1 t2

splinelen = lambda2-lambda1;

%% encode the constraints using rho-mu-tangent notation

rho=X(1);
mu=X(2);

t1=t1/norm(t1);
t2=t2/norm(t2);

a1 = t1(1)*rho;
a2 = (mu*t2(1)-a1)/(2*splinelen);

b1 = t1(2)*rho;
b2 = (mu*t2(2)-b1)/(2*splinelen);
c1 = t1(3)*rho;
c2 = (mu*t2(3)-c1)/(2*splinelen);


%% calculate the spline length
% 
X = lambda1:splinelen/100:lambda2;
Y = sqrt((a1+2*a2*X).^2 + (b1+2*b2*X).^2 + (c1+2*c2*X).^2);
calculated_splinelen = trapz(X,Y);


% calculated_splinelen = ((rho*(rho - mu*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3))))/(mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3))) + ...
%   (mu^2*(mu - rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3))))/(mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3))) + ...
%   (mu^2*rho^2*(t1(3)^2*(t2(1)^2 + t2(2)^2) - 2*t1(1)*t1(3)*t2(1)*t2(3) - 2*t1(2)*t2(2)*(t1(1)*t2(1) + t1(3)*t2(3)) + t1(2)^2*(t2(1)^2 + t2(3)^2) + ...
%      t1(1)^2*(t2(2)^2 + t2(3)^2))*log(mu^2 - mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)) + ...
%       mu*sqrt(mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)))))/...
%    (mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)))^(3/2) - ...
%   (mu^2*rho^2*(t1(3)^2*(t2(1)^2 + t2(2)^2) - 2*t1(1)*t1(3)*t2(1)*t2(3) - 2*t1(2)*t2(2)*(t1(1)*t2(1) + t1(3)*t2(3)) + t1(2)^2*(t2(1)^2 + t2(3)^2) + ...
%      t1(1)^2*(t2(2)^2 + t2(3)^2))*log(-rho^2 + mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)) + ...
%       rho*sqrt(mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)))))/...
%    (mu^2 + rho^2 - 2*mu*rho*(t1(1)*t2(1) + t1(2)*t2(2) + t1(3)*t2(3)))^(3/2))/2;


ceq = calculated_splinelen-splinelen;

c = [];
