function K=Kmax(X)

global Xp1 Xp2 t1 t2 lambda1 lambda2
% t11 = Xp1(1);
% t12 = Xp1(2);
% t13 = Xp1(3);
% t21 = Xp2(1);
% t22 = Xp2(2);
% t23 = Xp2(3);

splinelen = lambda2-lambda1;

%split up the input matrix into its constituents
rho=X(1);
mu=X(2);
t1=t1/norm(t1);
t2=t2/norm(t2);
%a0 = X(1,1);
a1 = t1(1)*rho;
a2 = (mu*t2(1)-a1)/(2*splinelen);
%b0 = X(2,1);
b1 = t1(2)*rho;
b2 = (mu*t2(2)-b1)/(2*splinelen);
%c0 = X(3,1);
c1 = t1(3)*rho;
c2 = (mu*t2(3)-c1)/(2*splinelen);


num = (2*(a2^2 + b2^2 + c2^2)^(3/2));
den = (a2^2*(b1^2 + c1^2) + (b2*c1 - b1*c2)^2 - 2* a1*a2*(b1*b2 + c1*c2) + a1^2*(b2^2 + c2^2));

%provide for the case of an identically straight line
EPSILON = 1e-8;
if num < EPSILON
    K = 0;
else
    K = num/den;
end

%use rho-mu-tangent parameterization and check that they produce the same
%result
% a1 -> rho*t11
% b1 -> rho*t12
% c1 -> rho*t13
% a2 -> (mu*t21 - rho*t11)/2
% b2 -> (mu*t22 - rho*t12)/2
% c2 -> (mu*t23 - rho*t13)/2
% 
% rho = a1/t11
% rho = b1/t12
% rho = c1/t13
% 
% mu = 2*a2+rho*t11/t21
% mu = 2*b2+rho*t12/t22
% mu = 2*c2+rho*t13/t23
% 
% 
% f2 = ((rho^2 * (t11^2 + t12^2 + t13^2) - 2*mu*rho*(t11*t21 + t12*t22 + t13*t23) + mu^2*(t21^2 + t22^2 + t23^2))^(3/2)) / (mu^2*rho^2*(t12^2*t21^2 + t13^2*t21^2 - 2*t11*t12*t21*t22 + t11^2*t22^2 + t13^2*t22^2 - 2*t11*t13*t21*t23 - 2*t12*t13*t22*t23 + t11^2*t23^2 + t12^2*t23^2))
% 
% diff = f-f2
% disp(diff);