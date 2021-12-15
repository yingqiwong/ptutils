function [X, Y, tri] = PlotTernaryBG (f, fName)
% 
% [X, Y, tri] = PlotTernaryBG (f, fName)
% 
% plots the ternary background for a three phase system
% 
% INPUTS
% f         phase fraction field [NPHS x N]
% fName     phase fraction names [NPHS x 1]
% 
% OUTPUTS
% X         x position of vertices [1 x N]
% Y         y position of vertices [1 x N]
% tri       triangle info [??]
% 
% YQW, 12 January 2021

% check size of f
if size(f,2)==3 && size(f,1)~=3
    f = f';
end

if nargin<2
    NPHS  = size(f,2);
    fName = strcat({'f'}, num2str((1:NPHS)', '%d'));
end

% load default plotting format
PlotDefFormat;

% prepare for plotting ternary plots
m      = 5;
n      = m-1;
grids  = linspace(0,1,m+1);
grids  = grids(1:end-1);
labels = num2str(grids(2:end)');

[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);
[X,Y]    =  terncoords(f(1,:),f(2,:),f(3,:));
tri      =  simpletri(size(f,1));

hold on;
% plot ternary axes
set(gca,'visible','off','defaultlinelinewidth',0.01);
maxz = 2;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;


% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,['$\phi^{' fName{1} '}$'],TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,['$\phi^{' fName{3} '}$'],TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,['$\phi^{' fName{2} '}$'],TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
colormap(ocean);
end