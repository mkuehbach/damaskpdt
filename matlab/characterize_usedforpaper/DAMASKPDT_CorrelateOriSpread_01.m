%% DAMASKPDT correlate orientation spread per mat point
%% using precomputed modal orientations per grain
% MTex v5.0.3, we utilize SI units
% Markus K\"uhbach, 2019/05/09 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;
    parm.plotresu = '-r600';
    parm.fontsz = 24;
    parm.fontnm = 'DejaVu Sans';
    parm.fontszax = 22;
    
    parm.opt_odf = 0;
	parm.opt_volfrac = 0;
    parm.opt_odfplot = 0; %1 - key section, 2 - classical 18, parula coloring
    

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
    
%% load user data
    parm.fns = dir(['*GrainOriSpatDistr*.bin']);
    parm.N = length(parm.fns);
    parm.NXYZ = int32(256^3);
    
%% run post-processing loop
for i=1:1:parm.N
    %fn = [parm.prefix '.Incr.' num2str(parm.i) '.GrainIDQuat' parm.suffix];
    parm.fn = parm.fns(i).name;
    display(parm.fn);
    parm.NR = str2double(extractBetween(parm.fn,'NR.','.NC'));
    parm.NC = str2double(extractBetween(parm.fn,'NC.','.bin'));
    parm.ID = str2double(extractBetween(parm.fn,'SimID.','.Rank'));
       
%% load raw heavy data from damaskpdt output
    display(['Loading ' parm.fn ' ...']);
    fid = fopen(parm.fn);
    GrainOriSpatDistr = fread(fid,[parm.NC, parm.NR],'double'); %implicit transpose
    fclose(fid);
    display([parm.fn ' loaded']);
    clearvars -except i cs ss parm GrainOriSpatDistr;

%% load precomputed modal orientations per grain
    parm.prefix = char(extractBetween(parm.fn,'DAMASKPDT.','.GrainOriSpatDistr'));
    %https://de.mathworks.com/help/matlab/matlab_prog/cell-arrays-of-strings.html
    parm.modalfn = ['DAMASKPDT.' parm.prefix '.GrainIDQuat.NR.' num2str(parm.NXYZ) '.NC.5.bin.mat'];
    % do not load acidentally also parm it would overwrite!
    load(parm.modalfn,'GrainIDsModalOris','cGrainIDsModalOris');
    
%% correlate disorientation to grain mean with distance
    %ecdf(GrainOriSpatDistr(2,:).*256)

    %walk through points and compute disorientation to mean
    GrainIDsDisoriToModal = zeros(4,length(GrainOriSpatDistr(1,:))); 
	% i) gid, 
	% ii) dist,
    % iii) disori with modal compute from damask passive no change previous and q_js also just taken from damask no change from passive to active pq^-1
	% iv) disori with modal computed after all damask passive to active then mtex mean ori and each q_j also taken from passive to active then angle( cmodal, cq_j) i.e. p^-1 q
    %assign grain Ids and distance
    GrainIDsDisoriToModal(1:2,:) = GrainOriSpatDistr(1:2,:);
	qin = GrainOriSpatDistr(3:6,:);
	clearvars GrainOriSpatDistr;

	% memory lean compromise of processing vectorized blocks of length nk
	%NK = 1000;
	%NN = length(GrainOriSpatDistr(1,:));
	%k = 1;
	%while k < NN
	%	if k + NK < NN 
	%		nk = NK;
	%	else
	%		nk = NN - k;
	%	end
	%	display(num2str(nk));
	%	
	%	ival = k:k+nk;
	%	p = 
	%	
	%	q = orientation(quaterion(qin(1,ival),qin(2,ival),qin(3,ival),qin(4,ival)),cs,ss); %original from damask passive
	%	cq = orientation(quaternion(qin(1,ival),-qin(2,ival),-qin(3,ival),-qin(4,ival)),cs,ss); % damask passive to active		
	%
	%end

    % build quaternion array measured
	q = orientation(quaternion(qin(1,:),qin(2,:),qin(3,:),qin(4,:)),cs,ss); %original from damask passive
	cq = orientation(quaternion(qin(1,:),-qin(2,:),-qin(3,:),-qin(4,:)),cs,ss); % damask passive to active		

    % build replicated array of modal orientations to speed this up
    gidmax = max(GrainIDsModalOris(1,:));
	cgidmax = max(cGrainIDsModalOris(1,:));
    % build first a O(1) time complex query array from which to pull precomputed ori class instances
    tmp = nan(4,1+gidmax);
	ctmp = nan(4,1+cgidmax);
    for j=1:1:length(GrainIDsModalOris(1,:))
        gid = GrainIDsModalOris(1,j);
        tmp(1:4,1+gid) = GrainIDsModalOris(2:5,j);
    end
	for j=1:1:length(cGrainIDsModalOris(1,:))
        cgid = cGrainIDsModalOris(1,j);
        ctmp(1:4,1+cgid) = cGrainIDsModalOris(2:5,j);
    end
    qtmp = orientation(quaternion(tmp),cs,ss);
    qctmp = orientation(quaternion(ctmp),cs,ss);
%    % in-memory in-place vectorized copying 
%    % orientation instances using constant time access querying field
    p = qtmp(1+GrainIDsDisoriToModal(1,:)); %modal computed with taking damask passive as is
    cp = qctmp(1+GrainIDsDisoriToModal(1,:)); %modal computed after conversion of damask passive to active
    clearvars gidmax cgidmx tmp ctmp j gid cgid qtmp qctmp;
    display('Disorientation preparation performed');
    
    % compute disorientation
    tic
    GrainIDsDisoriToModal(3,:) = angle(p,q) ./degree; %
	GrainIDsDisoriToModal(4,:) = angle(cp,cq) ./degree; %modal after active conversion, candidate from passive to active
    display('Matpoint disorientations to modal computed');
	toc
    clearvars ans cgidmax p cp q cq qin;
	%save([parm.fn '.mat'],'-v7.3');
 
 
    %ecdf(GrainIDsDisoriToModal(4,:))
    %hold on
    %ecdf(GrainIDsDisoriToModal(3,:))
    
    %% plot 2d histogram
    figure('Position',[100 100 1920 1920]);
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
    Xedges = [-Inf 0:0.1:24 Inf]; %distance in multiples of cells
    Yedges = [-Inf 0:0.5:63 Inf]; %disori
    h = histogram2(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(3,:), ...
        Xedges,Yedges, ...
        'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    colorbar;
    if contains(parm.fn,'.20')
        caxis([0.0 0.10])
    end
    if contains(parm.fn,'.100')
        caxis([0.0 0.020])
    end
    if contains(parm.fn,'.345')
        caxis([0.0 0.002])
    end
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
   
    pbaspect([1 1 0.5]) 
    axis tight
     xmin = 0;
     xmax = 24;
%     ymin = -300.0;
%     ymax = +50.0;
     xlim([xmin xmax]);
     xt = [xmin:(xmax-xmin)/12:xmax];
     xticks(xt);
%     ylim([ymin ymax]);
%     yt = [ymin:(ymax-ymin)/7:ymax];
%     yticks(yt);
%     %view(46,20);
    parm.fn = parm.fns(i).name;
    print(gcf,[parm.fn '.PasModalVsPas.png'],'-dpng',parm.plotresu); 
    close(gcf);
    

    %% plot 2d histogram
    figure('Position',[100 100 1920 1920]);
    box('on')
    %% variant 2 - build in histogram2d plot matlab
    Xedges = [-Inf 0:0.1:24 Inf]; %distance in multiples of cells
    Yedges = [-Inf 0:0.5:63 Inf]; %disori
    h = histogram2(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(4,:), ...
        Xedges,Yedges, ...
        'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    colorbar;
    if contains(parm.fn,'.20')
        caxis([0.0 0.10])
    end
    if contains(parm.fn,'.100')
        caxis([0.0 0.020])
    end
    if contains(parm.fn,'.345')
        caxis([0.0 0.002])
    end
    view(0,90);
    xlabel({'Distance to boundary'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    ylabel({'Disorientation to grain modal ori (°)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    set(gca,'FontSize',parm.fontsz,'FontName',parm.fontnm,'LineWidth',1.5)
    set(gcf,'PaperUnits','Inches')
    set(gcf,'PaperSize',[12/2.57 12/2.57])
    set(gcf,'color','w')
    pbaspect([1 1 0.5]) 
    axis tight
    xmin = 0;
    xmax = 24;
    xlim([xmin xmax]);
    xt = [xmin:(xmax-xmin)/12:xmax];
    xticks(xt);
    parm.fn = parm.fns(i).name;
    print(gcf,[parm.fn '.ActModalVsAct.png'],'-dpng',parm.plotresu); 
    close(gcf); 
    
    
    display([parm.fn ' processed completely']);
    
    
    
    clearvars -except i cs ss parm; 
end



















%% damaskpdt
% functionality to quantify the spatial distribution of state variable
% values with respect to the proximity at local interfaces
% Markus K�hbach, 2019/04/12 m.kuehbach at mpie.de
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

    
%% compute normalized mono-parameter eCDF for \rho for binning in distance classes with spacing equjent to the discretization
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