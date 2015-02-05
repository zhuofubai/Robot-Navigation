function [c,ceq] = nonlconst(X)
%constraints:
%starting point: X1
%ending point: Xp2
%starting tangent: Xp1
%ending tangent: Xp2
%spline length: splinelen

global splinelen X1 X2 Xp1 Xp2

ceq = [];
%add the tangent constraints

if ~isempty(Xp1)
    t11 = Xp1(1);
    t12 = Xp1(2);
    t13 = Xp1(3);
    ceq = [X(1,2) - t11;
           X(2,2) - t12;
           X(3,2) - t13];
end

if ~isempty(Xp2)
    t21 = Xp2(1);
    t22 = Xp2(2);
    t23 = Xp2(3);
    ceq = [ceq;
           X(1,2)+2*X(1,3)*splinelen - t21;
           X(2,2)+2*X(2,3)*splinelen - t22;
           X(3,2)+2*X(3,3)*splinelen - t23];
end

%add the starting and ending point constraints if they exist
if ~isempty(X1)
    x1 = X1(1);
    y1 = X1(2);
    z1 = X1(3);
    ceq = [ceq;
           X(1,1) - x1;
           X(2,1) - y1;
           X(3,1) - z1];
end

if ~isempty(X2)
    x2 = X2(1);
    y2 = X2(2);
    z2 = X2(3);
    ceq = [ceq;
           X(1,1) + X(1,2)*splinelen + X(1,3)*splinelen^2 - x2;
           X(2,1) + X(2,2)*splinelen + X(2,3)*splinelen^2 - y2;
           X(3,1) + X(3,2)*splinelen + X(3,3)*splinelen^2 - z2];
end

%find the length of the spline
c = abs((((X(1,2)*X(1,3) + X(2,2)*X(2,3) + X(3,2)*X(3,3))/(X(1,3)^2 + X(2,3)^2 + X(3,3)^2) + 2*splinelen)*sqrt(X(1,2)^2 + X(2,2)^2 + X(3,2)^2 + 4*X(1,2)*X(1,3)*splinelen + 4*X(2,2)*X(2,3)*splinelen + 4*X(3,2)*X(3,3)*splinelen + ...
     4*X(1,3)^2*splinelen^2 + 4*X(2,3)^2*splinelen^2 + 4*X(3,3)^2*splinelen^2) + ...
  ((X(1,3)^2*(X(2,2)^2 + X(3,2)^2) + (X(2,3)*X(3,2) - X(2,2)*X(3,3))^2 - 2*X(1,2)*X(1,3)*(X(2,2)*X(2,3) + X(3,2)*X(3,3)) + X(1,2)^2*(X(2,3)^2 + X(3,3)^2))*...
    log(2*(X(1,2)*X(1,3) + X(2,2)*X(2,3) + X(3,2)*X(3,3) + 2*X(1,3)^2*splinelen + 2*X(2,3)^2*splinelen + 2*X(3,3)^2*splinelen + sqrt(X(1,3)^2 + X(2,3)^2 + X(3,3)^2)*...
        sqrt(X(1,2)^2 + X(2,2)^2 + X(3,2)^2 + 4*X(1,2)*X(1,3)*splinelen + 4*X(2,2)*X(2,3)*splinelen + 4*X(3,2)*X(3,3)*splinelen + 4*X(1,3)^2*splinelen^2 + 4*X(2,3)^2*splinelen^2 + 4*X(3,3)^2*splinelen^2))))/...
   (X(1,3)^2 + X(2,3)^2 + X(3,3)^2)^(3/2))/4) - splinelen - 1e-4;%we only require the spline to be almost the right length