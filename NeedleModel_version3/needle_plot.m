function needle_plot(ABC,G)
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
        clf
        hold on
         plot3(G(:,3),G(:,1),G(:,2),'ro');
        plot3(ZData,XData,YData,'r');
        xlabel('Z-axis')
        ylabel('X-axis')
        zlabel('Y-axis')
        box on
        grid on