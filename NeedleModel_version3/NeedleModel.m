
clear all
close all
clf

%set up important variables
doVisualization = 1;
nControlPoints = 5;
nSteps = 20;
mu_delta = .1; %mean of the insertion distance distribution
sigma_delta = 0;%.01;
mu_d = 0; %mean of the normal distribution for the needle heading
sigma_d = 0;%.05; %std deviation


%set up an initially straight needle pointing straight along the z-axis
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

%declare some globals for the nonlinear constraint function to use
global lambda1 lambda2 t1 t2


fopts =  optimset('Display', 'off', 'Algorithm', 'sqp');%,'Display','Iter');%,'TolCon',0.); %options for fmincon



totalDelta = 0;


if doVisualization == 1
    %% plot
    XData = [];
    YData = [];
    ZData = [];
    for i=1:nControlPoints-1
        a0 = ABC(1,1,i);
        a1 = ABC(1,2,i);
        a2 = ABC(1,3,i);
        b0 = ABC(2,1,i);
        b1 = ABC(2,2,i);
        b2 = ABC(2,3,i);
        c0 = ABC(3,1,i);
        c1 = ABC(3,2,i);
        c2 = ABC(3,3,i);
        l = 0:.01:1;
        XData = [XData a0 + a1*l + a2*l.^2];
        YData = [YData b0 + b1*l + b2*l.^2];
        ZData = [ZData c0 + c1*l + c2*l.^2];
    end
    %plot each point and the spines connecting them
    clf
    hold on
    plot3(G(:,3),G(:,1),G(:,2),'bo');
    plot3(ZData,XData,YData,'b');
    xlabel('Z-axis')
    ylabel('X-axis')
    zlabel('Y-axis')
    box on
    grid on
    axis([0,nControlPoints+nSteps*mu_delta+.5,-1,1,-1,1]);
    title(sprintf('step %i, total insertion = %f',0,totalDelta));
    pause(.05);
end
for iStep=1:nSteps
    GLast = G;
    GpLast = Gp;
    ABCLast = ABC;
    
    delta = mu_delta + sigma_delta*randn; %choose delta from a normal distribution
    totalDelta = totalDelta+delta;
    
    %% find the spline from the old endpoint to the new endpoint
    
    cp = nControlPoints; %select the last control point as the current cp
    
    %first, find a new heading for the needletip
    deflection = mu_d + sigma_d * randn(1,3); %deviation of the new tangent from the old
%     deflection = [.05 0 0];
    
    Gpend = (Gp(cp,:) + deflection)/norm(Gp(cp,:) + deflection); %get the new unit-length tangent
    
    t1 = Gp(cp,:); %starting tangent, i.e., the tangent at the last control point
    t2 = Gpend; %ending tangent, i.e., the tangent at the new needletip
    lambda1 = 0;
    lambda2 = delta; %we're finding the delta-length spline from the endpoint to the new endpoint
    
    a0=G(cp,1);
    b0=G(cp,2);
    c0=G(cp,3);
    
    rhomuinit = [1 1];
    [rhomu,fval,exitflag] = fmincon(@Kmax,rhomuinit,[],[],[],[],[],[],@nonlconst,fopts);
    rho = rhomu(1);
    mu = rhomu(2);
    ABCend=[a0 rho*t1(1) (mu*t2(1)-rho*t1(1))/(2*(lambda2-lambda1));
        b0 rho*t1(2) (mu*t2(2)-rho*t1(2))/(2*(lambda2-lambda1));
        c0 rho*t1(3) (mu*t2(3)-rho*t1(3))/(2*(lambda2-lambda1))  ];
    
    
    %plug and chug to find the end point
    a0 = ABCend(1,1);
    a1 = ABCend(1,2);
    a2 = ABCend(1,3);
    b0 = ABCend(2,1);
    b1 = ABCend(2,2);
    b2 = ABCend(2,3);
    c0 = ABCend(3,1);
    c1 = ABCend(3,2);
    c2 = ABCend(3,3);
    lambda=delta;
    Gend = [a0 + a1*(lambda) + a2*(lambda)^2 b0 + b1*(lambda) + b2*(lambda)^2 c0 + c1*(lambda) + c2*(lambda)^2];
    
    %sanity checks
    t1_2 = [a1 b1 c1];
    t2_2 = [a1+2*a2*lambda b1+2*b2*lambda c1+2*c2*lambda];
    t1_2 = t1_2/norm(t1_2);
    t2_2 = t2_2/norm(t2_2);
    
    %now advance the last control point to the new end point
    G(cp,:) = Gend;
    Gp(cp,:) = Gpend/norm(Gpend);
    
    l = 0:.01:delta;
    XData = a0 + a1*l + a2*l.^2;
    YData = b0 + b1*l + b2*l.^2;
    ZData = c0 + c1*l + c2*l.^2;
    
    if doVisualization == 1
        %plot each point and the spines connecting them
        %     clf
        hold on
        plot3(ZData,XData,YData,'c');
        title('extension step');
        pause(.05);
    end
    
    %% now adjust the previous control points
    
    for nspline = nControlPoints-1:-1:1
        
        ABCcur = ABC(:,:,nspline);
        
        
        %advance the control point along the spline that starts at cp
        a0 = ABCcur(1,1);
        a1 = ABCcur(1,2);
        a2 = ABCcur(1,3);
        b0 = ABCcur(2,1);
        b1 = ABCcur(2,2);
        b2 = ABCcur(2,3);
        c0 = ABCcur(3,1);
        c1 = ABCcur(3,2);
        c2 = ABCcur(3,3);
        
        %get the control point at the base of the new spline
        X1 = [a0 + a1*delta + a2*delta^2 b0 + b1*delta + b2*delta^2 c0 + c1*delta + c2*delta^2];
        %get the control point at the end of the new spline
        X2 = G(nspline+1,:);
        %we know the endpoint tangent
        t2 = Gp(nspline+1,:);
        
        
        %now get t1 using klm parameterization
        k1 = X1(1);
        l1 = X1(2);
        m1 = X1(3);
        k2 = X2(1);
        l2 = X2(2);
        m2 = X2(3);
        
        xp2 = t2(1);
        yp2 = t2(2);
        zp2 = t2(3);
        
        k3 = -k1+k2-xp2;
        l3 = -l1+l2-yp2;
        m3 = -m1+m2-zp2;

        %solve for the tangents
        lambda=0;
        t1 = [-k1+k2+k3*(1-2*lambda) -l1+l2+l3*(1-2*lambda) -m1+m2+m3*(1-2*lambda)];
        t1 = t1/norm(t1);
        lambda = 1;
        t2 = [-k1+k2+k3*(1-2*lambda) -l1+l2+l3*(1-2*lambda) -m1+m2+m3*(1-2*lambda)]; %this should just be a sanity check
        t2 = t2/norm(t2);
        
        t2 = Gp(nspline+1,:); %we know this from the previous calulation, so it should match, and it does
        t2 = t2/norm(t2);
        
        a0 = X1(1);
        b0 = X1(2);
        c0 = X1(3);
        
        lambda1=0;
        lambda2=1;
        rhomuinit = [1 1];
        [rhomu,fval,exitflag] = fmincon(@Kmax,rhomuinit,[],[],[],[],[],[],@nonlconst,fopts);
        rho = rhomu(1);
        mu = rhomu(2);
        ABCnew=[a0 rho*t1(1) (mu*t2(1)-rho*t1(1))/(2*(lambda2-lambda1));
            b0 rho*t1(2) (mu*t2(2)-rho*t1(2))/(2*(lambda2-lambda1));
            c0 rho*t1(3) (mu*t2(3)-rho*t1(3))/(2*(lambda2-lambda1))  ];
        
        ABC(:,:,nspline) = ABCnew;
        G(nspline,:) = ABC(:,1,nspline);
        
        
        a0 = ABCnew(1,1);
        a1 = ABCnew(1,2);
        a2 = ABCnew(1,3);
        b0 = ABCnew(2,1);
        b1 = ABCnew(2,2);
        b2 = ABCnew(2,3);
        c0 = ABCnew(3,1);
        c1 = ABCnew(3,2);
        c2 = ABCnew(3,3);
        t1 = [a1+2*a2*delta b1+2*b2*delta c1+2*c2*delta];
        t1 = t1/norm(t1);
        if nspline > 1
            Gp(nspline,:) = t1;
        end
        
        %sanity checks
        lambda = 1;
        t1_2 = [a1 b1 c1];
        t2_2 = [a1+2*a2*lambda b1+2*b2*lambda c1+2*c2*lambda];
        t1_2 = t1_2/norm(t1_2);
        t2_2 = t2_2/norm(t2_2);
        lambda=0;
        X1_2 = [a0+a1*lambda+a2*lambda^2 b0+b1*lambda+b2*lambda^2 c0+c1*lambda+c2*lambda^2];
        lambda=1;
        X2_2 = [a0+a1*lambda+a2*lambda^2 b0+b1*lambda+b2*lambda^2 c0+c1*lambda+c2*lambda^2];
    end
    
    
    if doVisualization == 1
        %% plot
        XData = [];
        YData = [];
        ZData = [];
        for i=1:nControlPoints-1
            a0 = ABC(1,1,i);
            a1 = ABC(1,2,i);
            a2 = ABC(1,3,i);
            b0 = ABC(2,1,i);
            b1 = ABC(2,2,i);
            b2 = ABC(2,3,i);
            c0 = ABC(3,1,i);
            c1 = ABC(3,2,i);
            c2 = ABC(3,3,i);
            l = lambda1:.01:lambda2;
            XData = [XData a0 + a1*l + a2*l.^2];
            YData = [YData b0 + b1*l + b2*l.^2];
            ZData = [ZData c0 + c1*l + c2*l.^2];
        end
        %plot each point and the spines connecting them
        clf
        hold on
        plot3(Ginit(:,3),Ginit(:,1),Ginit(:,2),'bo');
        plot3(G(:,3),G(:,1),G(:,2),'ro');
        plot3(ZData,XData,YData,'r');
        xlabel('Z-axis')
        ylabel('X-axis')
        zlabel('Y-axis')
        box on
        grid on
        axis([0,nControlPoints+nSteps*mu_delta+.5,-1,1,-1,1]);
        title(sprintf('step %i, total insertion = %f',iStep,totalDelta));
        pause(.05);
    end
end
