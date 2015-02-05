clear all

St3(1).x=0.8708;St3(1).y=-2.501;St3(1).z=1.463;
St3(2).x=0.6228;St3(2).y=-2.1117;St3(2).z=1.673;
xsize=20;ysize=20;zsize=20;

[density A]=discrete_histgram(xsize,ysize,zsize, St3);
density
% figure,
% for i=1:1000
%     plot3(St2(i).x,St2(i).y,St2(i).z,'bo');
%     hold on
% end
% grid on