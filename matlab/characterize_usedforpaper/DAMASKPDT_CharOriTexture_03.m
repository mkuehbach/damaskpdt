%% DAMASKPDT characterize texture
% MTex v5.0.3, we utilize SI units
% Markus K\"uhbach, 2019/05/09 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;
    parm.plotresu = '-r600';
    parm.fontsz = 26;
    parm.fontnm = 'Calibri';
    parm.fontszax = 22;
    
%% load user data
    %global dbl;
    %global parms;
    parm.simid = 301402;
    parm.NXYZ = int32(256^3);
    
    parm.incr = '.Incr.345.1.345';
    parm.prefix = ['DAMASKPDT.SimID.' num2str(parm.simid) '.Rank.0' parm.incr];
    parm.suffix = '.NR.16777216.NC.5.bin';
    parm.opt_odf = 0;
    parm.opt_volfrac = 0;
    pamr.opt_odfplot = 0; %1 - key section, 2 - classical 18, parula coloring
 
    %parm.fn = [parm.prefix parm.suffix];
    %dbl.fns = dir([dbl.prefix '*.QuatCloud1.bin']); %'*.RA.12.RB.16.*.bin']);
    %dbl.n = length(dbl.fns);

%% run post-processing loop
for parm.i=345:1:345
    fn = [parm.prefix '.Incr.' num2str(parm.i) '.GrainIDQuat' parm.suffix];
    %display(fn);
    %parms.IncrID = str2double(extractBetween(fn,'Incr.','.Quat'));
    
%% load raw heavy data from damaskpdt output
    display(['Loading ' fn ' ...']);
    fid = fopen(fn);
    GrainIDsQuats = fread(fid,[5,parm.NXYZ],'double'); %implicit transpose
    fclose(fid);
    display([fn ' loaded']);
    clearvars -except parm GrainIDsQuats;

%% customize_mtex
    setMTEXpref('showCoordinates','on');
    setMTEXpref('FontSize',parm.fontsz);
    setMTEXpref('FontName',parm.fontnm);
    setMTEXpref('figSize','normal');
%% consistent_conventions
    % we redefine the MTex default coordinate system conventions from x2north
    % and z out of plane to x east and zinto plane which is the Setting 2 case
    % of TSL
    % https://github.com/mtex-toolbox/mtex/issues/56
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('zAxisDirection','intoPlane');
    setMTEXpref('defaultColorMap','parula');
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
    
    if parm.opt_odf == 1
%% build orientation list fast
        ori_in = orientation(quaternion(GrainIDsQuats(2:5,:)),cs,ss);
    
%% compute ODF
        parm.odf_resu = 10.0*degree;
        odf_apprx = calcODF(ori_in,'resolution',parm.odf_resu);
    
%% characterize ODF
        ['Texture index ' num2str(textureindex(odf_apprx))]
        parm.vol_scatter = 10.0*degree;
        %global volfrac;
        %volfrac.cube = volume(odf_apprx, orientation('Euler',0.0*degree,0.0*degree,0.0*degree,cs,ss), parms.vol_scatter)*100.0;
    
        %% write report file
        if parm.opt_volfrac == 1
            parm.fn_csv = [parm.prefix '.ODFVolFractions.csv'];
            fid = fopen(parm.fn_csv, 'wt');
            fprintf(fid, '%s\n',['Component;Scatter;VolumeFraction']);
            fprintf(fid, '%s\n',[';degree;percent']);
            fprintf(fid, '%s\n',['Component;Scatter;VolumeFraction']);
            fprintf(fid,'%s;%e;%e\n','Cube', parms.vol_scatter/degree, volfrac.cube );
            % report further
            fclose(fid);
        end
        
%% plot ODF classical or reduced, always triclinic
   
        if parm.opt_odfplot == 1
            % only key sections
            plot(odf_apprx,'phi2',[0.0 35.0 45.0]*degree,'contourf','silent');
            colorbar
            parm.fn_odf = [parm.mydir parm.prefix '.ODFDefResu.' num2str(parm.odf_resu/degree,'%.0f') '.KeyPhi2Sections.png'];
            parm.fn_odf = [parm.prefix '.ODFDefResu.' num2str(parm.odf_resu/degree,'%.0f') '.KeyPhi2Sections.png'];
        end
        if parm.opt_odfplot == 2
            % classical 5 degree sections
            plot(odf_apprx,'phi2','sections',18,'contourf','silent');
            colorbar
            parms.fn_odf = [parm.prefix '.ODFDefResu.' num2str(parm.odf_resu/degree,'%.0f') '.ClassicalPhi2.png'];
        end
       
        print(gcf,parm.fn_odf,'-dpng',parm.plotresu); 
        close(gcf);
    end

%% grain resolved orientation spread
    GrainIDsModalOris = unique(GrainIDsQuats(1,:));
    GrainIDsModalOris(2:5,:) = zeros(4,length(GrainIDsModalOris));
    
    %% extract all orientation values from matpoints of grain
    for j=1:1:length(GrainIDsModalOris(1,:))
        gid = GrainIDsModalOris(1,j);
        Quats = [];
        Quats = orientation(quaternion( ...
            GrainIDsQuats(2:5,GrainIDsQuats(1,:)==gid)),cs,ss);

        %% plotting options
        %    plotPDF(CurrOris,Miller(1,0,0,cs)) %PoleFigure 
        %    plotIPDF(CurrOris,vector3d(1,0,0)); %invPolefigure
        %    scatter(Quats,'rodrigues'); %in Rodrigues space

        %% mean of a list of orientations, principle axes and moments of inertia
        [m, lambda, V] = mean(Quats); %m - mean, lambda - principle moments of inertia (eigenvalues), V (principle axes of inertia ori, eigenvector), q2fz crystallographic equivalents projected to fundamental region
        GrainIDsModalOris(2:5,j) = [m.a;m.b;m.c;m.d];
        %% plot reference in there as well
            %plotPDF(m,Miller(1,0,0,cs))
            %hold on
            %scatter(m,'rodrigues','r')
        %% GROD and potentially add weights
        %% misorientation mis(o1,o2) = o1^-1 * o2;
        % average disorientation angle
%        val_km_vec = angle(Quats(1,:),m)./degree;
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
%         plot = 0;
%         if plot == 1
%             %modal phi1, Phi, phi2, lambda1234, mean disori angle against
%             %modal, q01, q25, q50, q75, q99
%             RESULTS(1+targetgid,1:3) = [m.phi1,m.Phi,m.phi2]./degree;
%             RESULTS(1+targetgid,4) = length(CurrOris(1,:));
%             RESULTS(1+targetgid,5:8) = lambda;
%             RESULTS(1+targetgid,9) = min(val_km_vec);
%             RESULTS(1+targetgid,10) = mean(val_km_vec); %mean orientation deviation
%             RESULTS(1+targetgid,11) = max(val_km_vec);
%             RESULTS(1+targetgid,12:16) = quantile(val_km_vec, [0.01, 0.25, 0.50, 0.75, 0.99]);
%         end 
 
        clearvars gid Quats m lambda V;
        display(num2str(j));
    end
    display('Grain modal orientations from matpoints computed');
    
%% correlate disorientation to grain mean with distance
    parm.suffix = '.NR.16777216.NC.6.bin';
    fn = [parm.prefix '.Incr.' num2str(parm.i) '.GrainOriSpatDistr' parm.suffix];
    display(['Loading ' fn ' ...']);
    fid = fopen(fn);
    GrainOriSpatDistr = fread(fid,[6,parm.NXYZ],'double'); %implicit transpose
    fclose(fid);
   
    %ecdf(GrainOriSpatDistr(2,:).*256)

    %walk through points and compute disorientation to mean
    GrainIDsDisoriToModal = zeros(3,length(GrainOriSpatDistr(1,:)));
    %assign grain Ids and distance
    GrainIDsDisoriToModal(1:2,:) = GrainOriSpatDistr(1:2,:);
    
    % build quaternion array measured
    qA = orientation(quaternion(GrainOriSpatDistr(3:6,:)),cs,ss);
    %qA_sw = ctranspose(qA); %inverse
    
    GrainOriSpatDistr(2:6,:) = [];
    
    % build replicated array of modal orientations to speed this up
    gidmax = max(GrainIDsModalOris(1,:));
    % build first a constant time querying structure from which to pull
    % precomputed orientation class instances
    tmp = nan(4,1+gidmax);
    for j=1:1:length(GrainIDsModalOris(1,:))
        gid = GrainIDsModalOris(1,j);
        tmp(1:4,1+gid) = GrainIDsModalOris(2:5,j);
    end
    qtmp = orientation(quaternion(tmp),cs,ss);
    % in-memory in-place vectorized copying 
    % orientation instances using constant time access querying field
    qB = qtmp(1+GrainOriSpatDistr(1,:));
    %qB_sw = ctranspose(qB);
    clearvars qtmp tmp gid j gidmax;
    display('Disorientation preparation performed');
    
    % compute disorientation
    tic
    GrainIDsDisoriToModal(3,:) = angle(qA,qB) ./degree;
    %GrainIDsDisoriToModal(4,:) = angle(qA_sw,qB_sw) ./degree;
    toc
    display('Matpoint disorientations to modal computed');
    
    %ecdf(GrainIDsDisoriToModal(4,:))
    %hold on
    %ecdf(GrainIDsDisoriToModal(3,:))
    
    %% plot 2d histogram
    figure('Position',[100 100 1920 1080]);
    %grid('on')
    box('on')
    %view(0,90)
    %view(58,24)

    %% plotting core
    %if opt.binning == 1
    %    %% variant 1 - 2d kernel density
    %    [bandwidth,density,X,Y]=kde2d(DAT);
    %    surf(X,Y,density,'LineStyle','none'), view([0,90])
    %    colormap parula, hold on, alpha(1.0)
    %    colorbar
    %end
    %if opt.binning == 2

    %% variant 2 - build in histogram2d plot matlab
    %h = histogram2(DAT(:,1),DAT(:,2),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');
    %h = histogram2(RHO(1,:)./parm.h,log10(RHO(2,:)),'XBinEdges',[1:1:15].*parm.h, 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    Xedges = [-Inf 0:0.1:20 Inf]; %distance in multiples of cells
    Yedges = [-Inf 0:0.5:63 Inf]; %disori
    h = histogram2(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(3,:), ...
        Xedges,Yedges, ...
        'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    colorbar;
    caxis([0.0 0.002])
    view(0,90);
    xlabel({'Distance to boundary'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    %if opt.dislotype == 1
    %ylabel({'\sigma_{33} (MPa)'},'FontSize',fontszcap,'FontName',fontnm);
    ylabel({'Disorientation to grain modal ori (°)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    %end
    %if opt.dislotype == 2
    %ylabel({'log_{10}(\Sigma\rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
    %end
    set(gca,'FontSize',parm.fontsz,'FontName',parm.fontnm,'LineWidth',1.5)
    set(gcf,'PaperUnits','Inches')
    set(gcf,'PaperSize',[12/2.57 12/2.57])
    set(gcf,'color','w')
    pbaspect([2 1 0.5])
%     xmin = 0.0;
%     xmax = 20;
%     ymin = -300.0;
%     ymax = +50.0;
%     xlim([xmin xmax]);
%     xt = [xmin:(xmax-xmin)/12:xmax];
%     xticks(xt);
%     ylim([ymin ymax]);
%     yt = [ymin:(ymax-ymin)/7:ymax];
%     yticks(yt);
%     %view(46,20);
    print(gcf,[fn '.png'],'-dpng',parm.plotresu); 
    close(gcf);
            
    %plot(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(3,:),'.')
      
end



















%% damaskpdt
% functionality to quantify the spatial distribution of state variable
% values with respect to the proximity at local interfaces
% Markus Kühbach, 2019/04/12 m.kuehbach at mpie.de
    clear;
    clc;
    digits 32;
    format long;
%% cases studied
    %SIMID = [0, 1, 2];
    SIMID = [40];
    %INCR = [0:10:100,114:14:240,255:15:390];
    INCR = [40];
    
%% user interaction
    parm.wkdir = ['D:/BIGMAX_RELATED/Paper/Paper13/bb_simulation/'];
    parm.NX = 75;
    parm.binning = 2;    %1 - ks2d, 2 - hist2d (binning no ksdensity!)
    %opt.dislotype = 1;  %1 - edges, 2 - dipoles, either case \sum all slipsystems per data point
    %if opt.dislotype == 1
    %    opt.dtype = 'EdgeDnsSpatDistr';
    %    parm.fns = dir([parm.prefix '*.' opt.dtype '*bin']);
    %else opt.dislotype == 2
    %    opt.dtype = 'DipoleDnsSpatDistr';
    %    parm.fns = dir([parm.prefix '*.' opt.dtype '*bin']);
    %end
        
    parm.h = 1.0/parm.NX;
    %parm.NC = 13; %distance + 12 primary slip system values
    parm.NC = 10; %distance + 3x3 tensor values

%% auto looping
    for c=1:1:length(SIMID(1,:))
        cc = SIMID(1,c);
        parm.prefix = [parm.wkdir 'DAMASKPDT.SimID.' num2str(cc)];
        for i=1:1:length(INCR(1,:))
            ii = INCR(1,i);
            fs = dir([parm.prefix '*.bin']);
            %find current file!
            for j=1:1:length(fs(:))
                if str2double(extractBetween(fs(j).name,'.Incr.','.Stress')) == ii
                    fn = fs(j).name;
                end
            end
            %fn = fs.name;
            %how many rows?
            NR = str2double(extractBetween(fn,'NR.', '.NC'));
            NC = parm.NC;
            
            %% open raw data
            fid = fopen([parm.wkdir fn]);
            DAT = fread(fid,[NC,NR],'double'); %implicit transpose
            fclose(fid);
            display(['Loaded ' fn]);
            %ecdf(DAT(1,:))
            %plot(DAT(1,:),DAT(10,:),'.')
            
            %histogram(DAT(1,:),20);
            %ecdf(DAT(1,:).*75.0)
            %hold on
            %ecdf(DAT(1,:).*75.0)
            %max(DAT(1,:))
            
            %% compute rho sum
            %RHO(1,:) = DAT(1,:); %distances to interface
            %RHO(2,:) = DAT(2,:);
            %for c=3:NC
            %    RHO(2,:) = RHO(2,:) + DAT(c,:);
            %end

            %% descriptive stuff

   
            clearvars -except SIMID INCR parm c cc i ii;
        end
    end
    
    'done'
    
    if opt.dislotype == 2
        print(gcf,[parms.fn_f64 '.png'],'-dpng','-r500'); close(gcf);
    end
%clearvars -except parm opt;

    
%% compute normalized mono-parameter eCDF for \rho for binning in distance classes with spacing equivalent to the discretization
    NX = parm.NXYZ^(1/3);
    SubStep = 1.0;
    clearvars refleft refright ref b binleft binright cand alpha h p SIGNIFICANCE;
    %% reference distribution
    refleft = 0.0/(SubStep*NX);
    refright = 1.0/(SubStep*NX);
    ref = log10(RHO(2, RHO(1,:) >= refleft & RHO(1,:) < refright));
    for b=1:1:(15*SubStep)
        clearvars alpha binleft binright cand h p;
        alpha = 0.05;
        binleft = double(b-1)/(SubStep*NX); %assuming implicit scaling
        binright = double(b)/(SubStep*NX);  
        
        cand = RHO(2, RHO(1,:) >= binleft & RHO(1,:) < binright);   
        display([num2str(b) '--> [' num2str(binleft) '  ' num2str(binright) ' ] ' num2str(length(cand(1,:))/parm.NXYZ)]);
        %display(median(log10(cand)))
    
        if ( ~isempty(cand) )
            %h --> 0 do not reject, h --> 1 reject null hypothesis
            [h, p] = kstest2(log10(cand),ref(1,:),'Alpha',alpha); %[h,p] if h == true rejecting null hypothesis at the 0.01 significance level
            SIGNIFICANCE(b,1) = parm.INCR(i);
            SIGNIFICANCE(b,2) = binleft;
            SIGNIFICANCE(b,3) = binright;
            SIGNIFICANCE(b,4) = double(h); %accept null hypothesis the same continuous distribution reject 1
            SIGNIFICANCE(b,5) = p; %sigl level
            SIGNIFICANCE(b,6) = double(alpha);
            SIGNIFICANCE(b,7) = length(cand);
            SIGNIFICANCE(b,8) = log10(min(cand));
            SIGNIFICANCE(b,9:15) = log10(quantile(cand,[0.01,0.1,0.25,0.5,0.75,0.9,0.99]));
            SIGNIFICANCE(b,16) = log10(max(cand));
        else
            break;
%             SIGNIFICANCE(b,1) = parm.INCR(i);
%             SIGNIFICANCE(b,2) = binleft;
%             SIGNIFICANCE(b,3) = binright;
%             SIGNIFICANCE(b,4) = double(1.0); %reject null hypothesis the same continuous distribution reject 1
%             SIGNIFICANCE(b,5) = NaN; %sigl level
%             SIGNIFICANCE(b,6) = double(alpha);
%             SIGNIFICANCE(b,7) = length(cand);
%             SIGNIFICANCE(b,8) = NaN;
%             SIGNIFICANCE(b,9:15) = NaN;
%             SIGNIFICANCE(b,16) = NaN;
        end
    end
    parms.fn_csv = [parms.fn_f64 '.Sgnfc.csv'];
    fid = fopen(parms.fn_csv, 'wt');
    fprintf(fid, '%s\n','IncrementID;BinLeft;BinRight;RejectNullHypo;PValue;SignificanceLevel;NValues;lg10RhoMin;lg10Rho1;lg10Rho10;lg10Rho25;lg10Rho50;lg10Rho75;lg10Rho90;lg10Rho99;lg10RhoMax');
    fclose(fid);
    dlmwrite(parms.fn_csv,SIGNIFICANCE,'-append');
    
%     % report further
%     formatspec = '%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e;%.3e\n';
%     for b=1:1:length(SIGNIFICANCE(:,1))
%         fprintf(fid, formatspec, ...
%             SIGNIFICANCE(b,1),  SIGNIFICANCE(b,2),  SIGNIFICANCE(b,3), ...
%              SIGNIFICANCE(b,4), SIGNIFICANCE(b,5), SIGNIFICANCE(b,6), ...
%               SIGNIFICANCE(b,7), SIGNIFICANCE(b,8), SIGNIFICANCE(b,9), ...
%                SIGNIFICANCE(b,10), SIGNIFICANCE(b,11), SIGNIFICANCE(b,12), ...
%                 SIGNIFICANCE(b,13), SIGNIFICANCE(b,14), SIGNIFICANCE(b,15), ...
%                  SIGNIFICANCE(b,16) );
%     end
%    fclose(fid);
end


% % %two-sided KolmogorovSmirnov reject against null hypothesis data are drawn
% % %from the same, non-parametric assumption free
% % 
% % %% plot the evolution of median and hypothesis testing results against strain
% % targetclass = 1;
% % % ##start building the figure
% % figure('Position',[100 100 1000 1000]);
% % %grid('on')
% % box('on')
% % %view(0,90)
% % fontszcap = 22;
% % fontszax = 22;
% % fontnm = 'Calibri';
% % for inc=1:length(SIGNIFICANCE(:))
% %     % different classes different color
% %     hold on
% %     plot(0.1,log10(SIGNIFICANCE(inc).Median(targetclass)),'.');
% %     %plot(SIGNIFICANCE(inc).Strain(targetclass),SIGNIFICANCE(inc).Median(targetclass),'.'); 
% % end
% % set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
% % set(gcf,'PaperUnits','Inches')
% % set(gcf,'PaperSize',[30 30])
% % set(gcf,'color','w')
% % pbaspect([2 2 2])
% % xmin = 0.0;
% % xmax = 0.32;
% % ymin = 15;
% % ymax = 16;
% % xlim([xmin xmax]);
% % xt = [xmin:(xmax-xmin)/8:xmax];
% % xticks(xt);
% % ylim([ymin ymax]);
% % yt = [ymin:(ymax-ymin)/10:ymax];
% % yticks(yt);
% % 
% % xlabel({'\epsilon_{vM}'},'FontSize',fontszcap,'FontName',fontnm);
% % if opt_dislotype == 1
% %     ylabel({'log_{10}(\Sigma_{s} \rho^s_{edge}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
% % end
% % if opt_dislotype == 2
% %     ylabel({'log_{10}(\Sigma_{s} \rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
% % end
% % 
% % if opt_dislotype == 1
% %     print(gcf,[fn '.EdgeTotal.png'],'-dpng','-r500'); close(gcf);
% % end
% % if opt_dislotype == 2
% %     print(gcf,[fn '.DipoleTotal.png'],'-dpng','-r500'); close(gcf);
% % end 