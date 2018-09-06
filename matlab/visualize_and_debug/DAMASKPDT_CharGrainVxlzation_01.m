%% Read AABB data and render 3d visualization of bounding boxes

% read DAMASKPDT.SimID.*.Incr.*.GrainVxlGrids.csv file information and
% generate 3d image of overlapping voxel container
clear;
clc;
digits 32;
format long;

SimID = 1005;
Incr = 0;
%[IDs,Boxes(:,1),Boxes(:,2),Boxes(:,3),Boxes(:,4),Boxes(:,5),Boxes(:,6)] =read_grainvxlgrids_csv('DAMASKPDT.SimID.0.Incr.0.GrainVxlGrids.csv');
fn_in = ['DAMASKPDT.SimID.' num2str(SimID) '.Incr.' num2str(Incr) '.GrainVxlGrids.csv'];
[IDs,Boxes] = read_grainvxlgrids_csv(fn_in);

% figure code

% global box
figure('Position',[100 100 2000 1000]);
aabb3d_plot([Boxes(1,1),Boxes(1,2),Boxes(1,3)],Boxes(1,4), ...
    Boxes(1,5),Boxes(1,6),'black',2.0,1.0); 

% local boxes
cmap = parula;
ngr = length(Boxes(:,1));
for gr=2:length(Boxes(:,1))
    cid = 1+int32(gr/ngr*63);
    aabb3d_plot([Boxes(gr,1),Boxes(gr,2),Boxes(gr,3)], ...
        Boxes(gr,4),Boxes(gr,5),Boxes(gr,6),cmap(cid,1:3),1.5,0.2); 
end
view(27,20);
% highlight particular box
gr = 3;
cid = 1+int32(gr/ngr*63);
aabb3d_plot([Boxes(gr,1),Boxes(gr,2),Boxes(gr,3)], ...
    Boxes(gr,4),Boxes(gr,5),Boxes(gr,6),cmap(cid,1:3),2.0,1.0);

%% modify axis appearance
fontszcap = 18;
fontnm = 'Calibri';
set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
set(gcf,'PaperUnits','Inches')
set(gcf,'PaperSize',[30 30])
set(gcf,'color','w')
box off;
grid off;
xlabel('X','FontSize',fontszcap);
ylabel('Y','FontSize',fontszcap);
zlabel('Z','FontSize',fontszcap);
print(gcf,[fn_in '.png'],'-dpng','-r600'); 
close(gcf);

%generate parula map for paraview via xml file
for i=1:length(cmap(:,1))
    display(['<Point x="' num2str(i) '" o="1" r="' num2str(cmap(i,1)), ...
        '" g="' num2str(cmap(i,2)) '" b="' num2str(cmap(i,3)) '"/>']);
end

%cube_plot([1,1,1],1,1,1,'r');
% Figure configurations
% Define the range of x-axis, y-axis, and z-axis in form of
% [xmin,xmax,ymin,ymax,zmin,zmax].
% axis([0,1,0,1,0,1]);
% Set the axis with equal unit.
%axis equal;
% Show grids on the plot


% Control the ticks on the axises
h = gca; % Get the handle of the figure
% h.XTick = 0:0.5:1;
% h.YTick = 0:0.5:1;
% h.ZTick = 0:0.5:1;
% Set the color as transparient
%material metal
%alpha('color');
%alphamap('rampup');
% Set the view point
%view(30,30);
%hold off;
% plot the figure in the form of eps with 600 ppi named 'filename'
% print(gcf,'-depsc2','-r600','filename.eps')