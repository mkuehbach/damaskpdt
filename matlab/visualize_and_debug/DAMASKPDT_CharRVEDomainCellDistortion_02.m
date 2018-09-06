%% plot the potentially triclinic bounding boxes deliminating the RVE datasets
% Markus K\"uhbach, 2018/05/16 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;

%% load grepped collection of RVE domain shape evolution data
    parms.simid = 10018;
    parms.INCR = int32([0,10,20,30,40,50,60,70,80,90,100,114,128,142,156,170,184,198,212,226,240,265,290,315,340,365,390,415]); %,440,465]);
    parms.prefix = ['E:/Paper13_OrientationGradients/cc_analysis/' num2str(parms.simid) '/'];
    parms.suffix = ['.RVEShapeEvo'];
    parms.shpfn =  [parms.prefix 'DAMASKPDT.SimID.' num2str(parms.simid) '.RVEShapeEvo.txt'];
    parms.base0fn = [parms.prefix 'DAMASKPDT.SimID.' num2str(parms.simid) '.RVEBase0Evo.txt'];
    parms.base1fn = [parms.prefix 'DAMASKPDT.SimID.' num2str(parms.simid) '.RVEBase1Evo.txt'];
    
   
%% load shape data, delete nan rows from grep format 
    SHP = read_greppedrveshape_01(parms.shpfn);
    SHP(any(isnan(SHP),2),:) = [];
  
    
%% load base vector data, delete nans
    BASE1 = read_greppedrveshape_01(parms.base1fn);
    BASE1 = BASE1(:,1:3);
    BASE1(any(isnan(BASE1),2),:) = [];

%% batch process strain increments
%figure('Position',[100 100 1000 1000]);
for incr=1:length(parms.INCR(1,:))
    clearvars -except parms SHP incr BASE1;
    sidxs = (incr-1)*4 +1;
    sidxe = sidxs+3;
    bidxs = (incr-1)*3 +1;
    bidxe = bidxs+2;
    
    DS = SHP(sidxs:sidxe,:); %zeros(4,24);
    BS = BASE1(bidxs:bidxe,:);
    %NAMES = strings(4,1);
    %figure('Position',[100 100 1000 1000]);
    cmap = parula;
    
    for i=1:4
        % there is still an ordering problem within DAMASKPDT different for each
        % box and incorrect h,k,l for the

        ver = [     DS(i,1) DS(i,2), DS(i,3);
                    DS(i,4) DS(i,5), DS(i,6);
                    DS(i,7) DS(i,8), DS(i,9);
                    DS(i,10) DS(i,11), DS(i,12);
                    DS(i,13) DS(i,14), DS(i,15);
                    DS(i,16) DS(i,17), DS(i,18);
                    DS(i,19) DS(i,20), DS(i,21);
                    DS(i,22) DS(i,23), DS(i,24)  ];

        %  Define the faces of the unit cubic
        fac = [1 2 3 4;
            4 3 5 6;
            6 7 8 5;
            1 2 8 7;
            6 7 1 4;
            2 3 5 8];

        bx = [ver(:,1),ver(:,2),ver(:,3)];
        %plot in existent one  
        hold on;
        linethickness = 2.0;
        linealpha = 1.0;
        if i == 1
           edgecolor = 'black';
           edgecolor = cmap(1,:);
        end
        if i == 2
          edgecolor = 'red';
          edgecolor = cmap(11,:);
        end
        if i == 3
          edgecolor = 'blue';
          edgecolor = cmap(21,:);
        end
        if i == 4
            edgecolor = 'green';
            edgecolor = cmap(31,:);
        end

        %if i == 1
        %    hold on
            patch('Faces',fac,'Vertices',bx,'EdgeColor',edgecolor,'FaceAlpha',0.0,'LineWidth',linethickness,'EdgeAlpha',linealpha); %transparent faces
        %end
    end

%end
    %% add base vector
    %https://de.mathworks.com/matlabcentral/fileexchange/3345-plot-arrowhead
    % input: x1,y1 - starting point
    % x2,y2 - end point
    % options - come as pairs of "property","value" as defined for "line" and "patch"
    % controls, see matlab help for listing of these properties.
    % note that not all properties where added, one might add them at the end of this file.
    %
    % additional options are:
    % 'headwidth': relative to complete arrow size, default value is 0.07
    % 'headheight': relative to complete arrow size, default value is 0.15
    % (encoded are maximal values if pixels, for the case that the arrow is very long)
    %
    % output: handles - handles of the graphical elements building the arrow
    %
    % Example: plot_arrow( -1,-1,15,12,'linewidth',2,'color',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5] );
  
    %plot_arrow( 0,0,1.0,2.0,'linewidth',2,'headwidth',0.25,'headheight',0.33 );
    
    % plot 3 base vector deformed config
    for j=1:3
        mArrow3( [0.0 0.0 0.0], [BS(1,j) BS(2,j) BS(3,j)],'facealpha', 1.0, 'color', 'black', 'stemWidth', 0.01); 
    end
    % plot_arrow; % will launch demo

    fn_in = ['DAMASKPDT.SimID.' num2str(parms.simid) '.Incr.' num2str(parms.INCR(incr))];
    view(27,20);
    fontszcap = 20;
    fontnm = 'Calibri';
    set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
    set(gcf,'PaperUnits','Inches')
    set(gcf,'PaperSize',[30 30])
    set(gcf,'color','w')
    box('on'); %ff;
    grid off;
    xlabel('X','FontSize',fontszcap);
    ylabel('Y','FontSize',fontszcap);
    zlabel('Z','FontSize',fontszcap);
    %% ##MK:: use that box DS(4,:) has maximum values
    xmin = -1.5;
    xmax = +2.5;
    ymin = -1.5;
    ymax = +2.5;
    zmin = -1.5;
    zmax = +2.5;
    xlim([xmin xmax]);
    xt = [xmin:(xmax-xmin)/8:xmax];
    ylim([ymin ymax]);
    yt = [ymin:(ymax-ymin)/8:ymax];
    zlim([zmin zmax]);
    zt = [zmin:(zmax-zmin)/8:zmax];

    print(gcf,[parms.prefix fn_in parms.suffix '.png'],'-dpng','-r600'); 
    close(gcf);
    display(fn_in);
end
