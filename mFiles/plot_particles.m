clc
clear all

restoredefaultpath
addpath('H:\PhD\simulation_model\MC_interact\src\build','-end');

coords = readmatrix('coords.txt');
data = readmatrix('data.txt');

rHydr = data(:,1);
rMag = data(:,2);
lBox = data(1,3);
numbPart = data(2,3);
% rNebr = data(3,3);
% rCut = data(4,3);
% lCell = data(5,3);
% nCells = data(6,3);
L = length(coords);

pos = coords(:,1:3);
posMm = coords(:,7:9);
posEa = coords(:,10:12);

time = coords(1:numbPart:end,13);

fig = figure(1);
[X,Y,Z] = sphere(30); %coordinates of sphere

j = L/numbPart;

for i = 0:10:j-1
    clf(fig)
    idx = i*numbPart;
    
    for ip = 1:numbPart
        xHydr = X.*rHydr(ip);
        yHydr = Y.*rHydr(ip);
        zHydr = Z.*rHydr(ip);
        
        sHydr = surf(xHydr + pos(ip + idx,1),yHydr + pos(ip + idx,2),zHydr + pos(ip + idx,3),'FaceColor','c','FaceAlpha',1);
        hold on
        sHydr.EdgeColor = 'none';        
    end
    
    light
    %plot arrows for easy axis and magnetic moment
    ran = (1 + i*numbPart):(numbPart*(i+1));
    qEa = quiver3(pos(ran,1),pos(ran,2),pos(ran,3),posEa(ran,1),posEa(ran,2),posEa(ran,3),'k'); %arrows of easy axes
    qEa.LineWidth = 1;
    qEa.AutoScaleFactor = 0.7;
    
    qMm = quiver3(pos(ran,1),pos(ran,2),pos(ran,3),posMm(ran,1),posMm(ran,2),posMm(ran,3),'m'); %arrows of moment
    qMm.LineWidth = 1;
    qMm.AutoScaleFactor = 0.7;
    
    a = lBox/2;
    
    x = [-a -a];
    y = [-a a];
    z = [-a -a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a -a];
    y = [-a a];
    z = [a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [a a];
    y = [-a a];
    z = [a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [a a];
    y = [-a a];
    z = [-a -a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a a];
    y = [a a];
    z = [-a -a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a a];
    y = [a a];
    z = [a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [a a];
    y = [a a];
    z = [-a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a -a];
    y = [a a];
    z = [-a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a a];
    y = [-a -a];
    z = [-a -a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a a];
    y = [-a -a];
    z = [a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [a a];
    y = [-a -a];
    z = [-a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    x = [-a -a];
    y = [-a -a];
    z = [-a a];
    line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
    
    fig.Position(1:4) = [300 70 800 700];
    ax = gca;
    ax.PlotBoxAspectRatio = [1 1 1];
    axis([-a a -a a -a a])
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    
    xlabel('x-position in m')
    ylabel('y-position in m')
    zlabel('z-position in m')
    
   title("simulation progress " + time(i+1) + "s");
    
%     if i == 0
%         x = linspace(-a,a,2);
%         for cz = 0:nCells
%             z = ones(1,2)*cz*lCell - a;
%             y = ones(1,2)*lCell*nCells - 3*a;
%             line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
%         end
%         z = linspace(-a,a,2);
%         for cx = 0:nCells
%             x = ones(1,2)*cx*lCell - a;
%             y = ones(1,2)*lCell*nCells - 3*a;
%             line(x,y,z,'Color','#D95319','LineStyle','--','LineWidth',1)
%         end
%         
%         %plot circles
%         theta = linspace(0,2*pi,50);
%         xCircleCut = rCut*cos(theta);
%         zCircleCut = rCut*sin(theta);
%         xCircleShell = rNebr*cos(theta);
%         zCircleShell = rNebr*sin(theta);
%         yCircle = zeros(1,50) - a;
%         plot3(xCircleCut,yCircle,zCircleCut,'LineStyle','-','Color','#D95319')
%         plot3(xCircleShell,yCircle,zCircleShell,'LineStyle','-','Color','#D95319')
%     end
    
    drawnow;
    pause(0.1);
end
