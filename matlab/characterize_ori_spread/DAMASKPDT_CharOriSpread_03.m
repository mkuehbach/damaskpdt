%% characterize orientation spread supplemental to damaskpdt via MTex v4.5.2, we utilize SI units
% Markus K\"uhbach, 2018/05/03 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;

%% load user data
    %global dbl;
    %global parms;
    dbl.simid = 10136;
    dbl.INCR = int32([0,10,20,30,40,50,60,70,80,90,100,114,128,142,156,170,184,198,212,226,240,265,290,315,340,365,390,415]); %,440,465]);
    %dbl.simid = 666;
    %dbl.INCR = int32([415]);
    dbl.prefix = ['E:/Paper13_OrientationGradients/cc_analysis/1000/' num2str(dbl.simid) '/'];  %'/DAMASKPDT.SimID.' num2str(dbl.simid) '.Incr.0.1.415.Incr.'];
    %dbl.prefix = 'Z:/damaskpdt/run/500g128x128x128_compZ/';
    dbl.fns = dir([dbl.prefix '*.OriSpatDistr*.bin']);
    dbl.n = length(dbl.fns);
%% run post-processing loop
    plot = 1; 
    
for i=1:2 %1:length(dbl.INCR(1,:)) % maximum 415 -- 23:1:23 %dbl.n:-1:1
    %%%fn = dbl.fns(i).name;
    %%%parms.IncrID = dbl.INCR(i); %str2double(extractBetween(fn,'Incr.','.Quat'));
       
%% input file details
    %dbl.simid = 0;
    %parms.IncrID = 0;
    parms.NXYZ = 128^3;
    %%%parms.mydir = dbl.prefix; %'X:\damaskpdt\run\500g128x128x128_compZ\';
    %%%parms.prefix  = ['DAMASKPDT.SimID.' num2str(dbl.simid) '.Incr.0.1.415.Incr.' num2str(parms.IncrID)];
    %parms.fn_u32 = [parms.mydir parms.prefix '.' 'QuatCloudGID.bin'];
    %parms.fn_f32 = [parms.mydir parms.prefix '.' 'QuatCloudORI.bin']; 
    parms.fn_f64 = [dbl.prefix dbl.fns(i).name];
    
%% load raw heavy data from damaskpdt output
    display(['Loading ' parms.fn_f64 ' ...']);
    fid = fopen(parms.fn_f64);
    QUAT = fread(fid,[6,parms.NXYZ],'double'); %implicit transpose
    fclose(fid);
    
    %Distances of integration points row 1
    %GrainIDs row 2
    %rest components of active quaternion ori 0,1,2,3
      
    clearvars -except dbl parms i QUAT plot;

%% customize_mtex
    parms.fontsz = 26;
    setMTEXpref('showCoordinates','on');
    setMTEXpref('FontSize',parms.fontsz);
    setMTEXpref('figSize','normal');
%% consistent_conventions
    % we redefine the MTex default coordinate system conventions from x2north
    % and z out of plane to x east and zinto plane which is the Setting 2 case
    % of TSL
    % https://github.com/mtex-toolbox/mtex/issues/56
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('zAxisDirection','intoPlane');
    %by setting z such, we interpret and define from now on a handedness --- a right-handed --- namely
    %by extending the 2D geometrical SEM/EBSD csys into consistency with the Bunge reference
    getMTEXpref('EulerAngleConvention')
    getMTEXpref('xAxisDirection') %check how MTex is setup now
    getMTEXpref('zAxisDirection')
    %right-handed coordinate system
    
%% symmetries definition
    cs = crystalSymmetry('cubic');
    ss = specimenSymmetry('triclinic');
    display('MTex initialized and fcc -1 symmetries loaded...');

    gidmax = max(QUAT(2,:));
    %% extract all quaternions for all ips of targetgrain
    
    %% finally build relation between individual ip
    % distance and local disorientation angle to grain modal ori
    THETA = QUAT(1:3,:); %##MK::dummy copy over q0 will later be overwritten by correct value!
    THETA(3,:) = -inf;
    
    for targetgid=0:1:gidmax
        CurrCloud = QUAT(1:6,QUAT(2,:) == targetgid);
        
        % assure that we can copy order preserving and vectorized
        %here = find(THETA(2,:) == targetgid);
        %WTF = CurrCloud(1,:);
        %for j=1:length(WTF(1,:))
        %    WTF(1,j) = j*1000;
        %end
        %THETA(3,THETA(2,:) == targetgid) = WTF;
     
        % build orientation list fast
        CurrOris = orientation(quaternion( ...
            CurrCloud(3,:),CurrCloud(4,:),CurrCloud(5,:),CurrCloud(6,:)),cs,ss);

        %% plotting options
        %    plotPDF(CurrOris,Miller(1,0,0,cs)) %PoleFigure 
        %    plotIPDF(CurrOris,vector3d(1,0,0)); %invPolefigure
        %    scatter(CurrOris,'rodrigues'); %in Rodrigues space

        %% mean of a list of orientations, principle axes and moments of inertia
            [m, q2fz, lambda, V] = mean(CurrOris); %m - mean, lambda - principle moments of inertia (eigenvalues), V (principle axes of inertia ori, eigenvector), q2fz crystallographic equivalents projected to fundamental region

        %% plot reference in there as well
            %plotPDF(m,Miller(1,0,0,cs))
        %% GROD and potentially add weights
            %#########
        %% misorientation mis(o1,o2) = o1^-1 * o2;
        % average disorientation angle
        val_km_vec = angle(CurrOris(1,:),m)./degree;
        %val_mk_vec = angle(m, CurrOris(1,:)) ./degree;
        %for k=1:length(CurrOris(1,:))
        %    val_km_seq(1,k) = angle(CurrOris(1,k),m) / degree;
        %    val_mk_seq(1,k) = angle(m,CurrOris(1,k)) / degree;
        %end
        %compare
        %sum(abs(val_km_vec-val_km_seq))
        %mean(val_km_vec)

        %#########
        %q1 = CurrOris(1);
        %q2 = orientation.brass(cs,ss);
        %ans = inv(q1) * q2;
        
        if plot == 1
            %modal phi1, Phi, phi2, lambda1234, mean disori angle against
            %modal, q01, q25, q50, q75, q99
            RESULTS(1+targetgid,1:3) = [m.phi1,m.Phi,m.phi2]./degree;
            RESULTS(1+targetgid,4) = length(CurrOris(1,:));
            RESULTS(1+targetgid,5:8) = lambda;
            RESULTS(1+targetgid,9) = min(val_km_vec);
            RESULTS(1+targetgid,10) = mean(val_km_vec); %mean orientation deviation
            RESULTS(1+targetgid,11) = max(val_km_vec);
            RESULTS(1+targetgid,12:16) = quantile(val_km_vec, [0.01, 0.25, 0.50, 0.75, 0.99]);
        end 
        
        %COLLECT(1+targetgid,i) = mean(val_km_vec);
        
        %copy over order preserving disori angle to modal for
        %characterizing spread
        THETA(3,THETA(2,:) == targetgid) = val_km_vec;
     
        clearvars CurrCloud CurrOris m q2fz lambda V val_km_vec;
        display(num2str(targetgid));
    end

    %% report results
    if plot == 1
        parms.fn_csv = [parms.fn_f64 '.OriSpread.csv'];
        fid = fopen(parms.fn_csv, 'wt');
        fprintf(fid, '%s\n', ['ModalPhi1;ModalPhi;ModalPhi2;ValueCount;L4;L3;L2;L1;MinDisoriAngle;MeanDisoriAngleToModal;MaxDisoriAngle;Quantile1;Quantile25;Quantile50;Quantile75;Quantile99']);
        fprintf(fid, '%s\n',['degree;degree;degree;1;1;1;1;1;degree;degree;degree;degree;degree;degree']);
        fprintf(fid, '%s\n', ['ModalPhi1;ModalPhi;ModalPhi2;ValueCount;L4;L3;L2;L1;MinDisoriAngle;MeanDisoriAngleToModal;MaxDisoriAngle;Quantile1;Quantile25;Quantile50;Quantile75;Quantile99']);
        % report further
        formatspec = '%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e\n';
        for j=0:1:gidmax
            fprintf(fid, formatspec ,RESULTS(1+j,:));
        end
        fclose(fid);
        %cdfplot(RESULTS(:,8))
    end
    
    
    if plot == 1
        %% visualize orientation spread
        figure('Position',[100 100 2000 1000]);
        %grid('on')
        box('on')
        %view(0,90)
        fontszcap = 22;
        fontszax = 22;
        fontnm = 'Calibri';
        %build in histogram2d plot matlab
        pxh = 1.0 / (parms.NXYZ^(1/3));
        h = histogram2(THETA(1,:)/pxh,THETA(3,:),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
        colorbar;

        xlabel({'Distance (\cdot h)'},'FontSize',fontszcap,'FontName',fontnm);
        ylabel({'\Theta (°)'},'FontSize',fontszcap,'FontName',fontnm);

        set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
        set(gcf,'PaperUnits','Inches')
        set(gcf,'PaperSize',[30 30])
        set(gcf,'color','w')
        pbaspect([1 2 0.5])


        xmin = 0.0; %resolution in effect of grid size
        xmax = max(THETA(1,:)) / pxh + 1;
        ymin = 0.0; %disorientation
        ymax = 62.8;
        xlim([xmin xmax]);
        xt = [xmin:(xmax-xmin)/8:xmax];
        xticks(xt);
        ylim([ymin ymax]);
        yt = [ymin:(ymax-ymin)/8:ymax];
        yticks(yt);
        %ytickformat('$%,.1f');
        xtickformat('%.0f');
        ytickformat('%.1f');
        view(90,90);

        print(gcf,[parms.fn_f64 '.Hist2.png'],'-dpng','-r600'); close(gcf);
    end

    clearvars -except dbl i COLLECT plot;
    display(['Increment ' num2str(i) ' processed']);
end


%% collect values from OriSpread to plot evolution as a function of strain
%% this does not allow to place the box plots at the correspondingly numerically
% intepreted position on an x-axis showing true strain, one really has to
% report rowwise strains and corresponding quantile values
% these.fns = dir([dbl.prefix '*.OriSpread.csv']);
% clearvars COLLECT;
% for f=1:1:length(these.fns(:,1))
%     CurrIncr = str2double(extractBetween(these.fns(f).name,'415.Incr.','.OriSpat'));
%     DAT = read_orispread_mtex([dbl.prefix these.fns(f).name]);
%     %COLLECT(:,find(dbl.INCR == CurrIncr)) = DAT(:,10);
% end
% nl = length(COLLECT(:,1));
% MATRIX(2:1+nl,:) = COLLECT;
% % strain values from flowcurve
% flowcurve.fns = dir([dbl.prefix '*Flowcurve.csv']);
% [INCR, EPSVM, SIGMAVM] = read_flowcurve_damaskpdt([dbl.prefix flowcurve.fns(1).name]);
% MATRIX(1,:) = reshape(EPSVM,[1,length(EPSVM(:,1))]);
% csvwrite([dbl.prefix 'MeanGrainDisoriAngleVsTrueStrain.csv'],MATRIX);

%% but this works as intended
    these.fns = dir([dbl.prefix '*.OriSpread.csv']);
    % strain values from flowcurve
    flowcurve.fns = dir([dbl.prefix '*Flowcurve.csv']);
    [COLLECT(:,1), COLLECT(:,2), COLLECT(:,3)] = read_flowcurve_damaskpdt([dbl.prefix flowcurve.fns(1).name]);
    for f=1:1:length(these.fns(:,1))
        CurrIncr = str2double(extractBetween(these.fns(f).name,'415.Incr.','.OriSpat'));
        DAT = read_orispread_mtex([dbl.prefix these.fns(f).name]);
        COLLECT(find(dbl.INCR == CurrIncr),4:10) = quantile(DAT(:,10),[0.01,0.10,0.25,0.5,0.75,0.9,0.99]);
    end
    csvwrite([dbl.prefix 'MeanGrainDisoriAngleVsTrueStrain.csv'],COLLECT);

