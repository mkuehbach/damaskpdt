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
    parm.NX = 256;
    parm.NXYZ = parm.NX^3;
    
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
    
%% compute disorientation using MTex
    tic
    GrainIDsDisoriToModal(3,:) = angle(p,q) ./degree; %
	GrainIDsDisoriToModal(4,:) = angle(cp,cq) ./degree; %modal after active conversion, candidate from passive to active
    display('Matpoint disorientations to modal computed');
	toc
    clearvars ans cgidmax p cp q cq qin;
	
%% quantify quantiles of distance-dependent sub-distributions
    nquants = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95,0.99];
    nq = length(nquants);
    dbins = [0.0:0.2:24.0]./parm.NX; %end values
    Nbins = length(dbins)-1;
        
    %% build reference for KS Test
    %dright = 0.0/parm.NX;
    %dleft =  1.0/parm.NX;
    %refcand = [];
    %refvalues = [];
    %refcand = Distances(1,:) > dright & Distances(1,:) <= dleft; 
    %refvalues = StressSpatDistr(c,refcand);
   
%% PasModalVsPas
    QUANTMATRIX = nan(Nbins,1+1+1+nq+1); %+1+1+1+1); %binid, cnts, min, quant values, max
    %accept yes/no 0.05, p@0.05, accept yes/no @0.01, p@0.01
    %figure
    for d=1:1:Nbins
		clearvars cand values;
		dright = dbins(d);
		dleft = dbins(d+1);
		QUANTMATRIX(d,1) = dleft * parm.NX;
		cand = GrainIDsDisoriToModal(2,:) > dright & GrainIDsDisoriToModal(2,:) <= dleft; 
		QUANTMATRIX(d,2) = sum(cand);
		values = GrainIDsDisoriToModal(3,cand);
		if ~isempty(values)
			QUANTMATRIX(d,2+1) = min(values);
			QUANTMATRIX(d,2+1+1:2+1+nq) = quantile(values,nquants);
			QUANTMATRIX(d,2+1+nq+1) = max(values);
		end
		display(d);
	end	
    parm.csv = [parm.fn '.PasModalVsPas.csv'];
	fid = fopen(parm.csv, 'wt');
    fprintf(fid, '%s\n','DistanceBinEndPixel;Counts;Min;0.01;0.05;0.25;0.5;0.75;0.95;0.99;Max'); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
    for d=1:1:Nbins
		fprintf(fid, '%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f\n',QUANTMATRIX(d,:)); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
    end
    fclose(fid);
	
%% ActModalVsAct
	QUANTMATRIX = nan(Nbins,1+1+1+nq+1); %+1+1+1+1); %binid, cnts, min, quant values, max
    %accept yes/no 0.05, p@0.05, accept yes/no @0.01, p@0.01
    %figure
    for d=1:1:Nbins
		clearvars cand values;
		dright = dbins(d);
		dleft = dbins(d+1);
		QUANTMATRIX(d,1) = dleft * parm.NX;
		cand = GrainIDsDisoriToModal(2,:) > dright & GrainIDsDisoriToModal(2,:) <= dleft; 
		QUANTMATRIX(d,2) = sum(cand);
		values = GrainIDsDisoriToModal(4,cand);
		if ~isempty(values)
			QUANTMATRIX(d,2+1) = min(values);
			QUANTMATRIX(d,2+1+1:2+1+nq) = quantile(values,nquants);
			QUANTMATRIX(d,2+1+nq+1) = max(values);
		end
		display(d);
	end	
	parm.csv = [parm.fn '.ActModalVsAct.csv'];
	fid = fopen(parm.csv, 'wt');
    fprintf(fid, '%s\n','DistanceBinEndPixel;Counts;Min;0.01;0.05;0.25;0.5;0.75;0.95;0.99;Max'); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
    for d=1:1:Nbins
		fprintf(fid, '%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f\n',QUANTMATRIX(d,:)); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
    end
    fclose(fid);
    
    display([parm.fn ' processed completely']);
    clearvars -except i cs ss parm;
end
	








	%save([parm.fn '.mat'],'-v7.3');
 
 
    %ecdf(GrainIDsDisoriToModal(4,:))
    %hold on
    %ecdf(GrainIDsDisoriToModal(3,:))
    
    % %% plot 2d histogram
    % figure('Position',[100 100 1920 1920]);
    % %grid('on')
    % box('on')
    % %view(0,90)
    % %view(58,24)

    % %% plotting core
    % %if opt.binning == 1
    % %    %% variant 1 - 2d kernel density
    % %    [bandwidth,density,X,Y]=kde2d(DAT);
    % %    surf(X,Y,density,'LineStyle','none'), view([0,90])
    % %    colormap parula, hold on, alpha(1.0)
    % %    colorbar
    % %end
    % %if opt.binning == 2

    % %% variant 2 - build in histogram2d plot matlab
    % %h = histogram2(DAT(:,1),DAT(:,2),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');
    % %h = histogram2(RHO(1,:)./parm.h,log10(RHO(2,:)),'XBinEdges',[1:1:15].*parm.h, 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    % Xedges = [-Inf 0:0.1:24 Inf]; %distance in multiples of cells
    % Yedges = [-Inf 0:0.5:63 Inf]; %disori
    % h = histogram2(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(3,:), ...
        % Xedges,Yedges, ...
        % 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    % colorbar;
    % if contains(parm.fn,'.20')
        % caxis([0.0 0.10])
    % end
    % if contains(parm.fn,'.100')
        % caxis([0.0 0.020])
    % end
    % if contains(parm.fn,'.345')
        % caxis([0.0 0.002])
    % end
    % view(0,90);
    % xlabel({'Distance to boundary'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    % %if opt.dislotype == 1
    % %ylabel({'\sigma_{33} (MPa)'},'FontSize',fontszcap,'FontName',fontnm);
    % ylabel({'Disorientation to grain modal ori (°)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    % %end
    % %if opt.dislotype == 2
    % %ylabel({'log_{10}(\Sigma\rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
    % %end
    % set(gca,'FontSize',parm.fontsz,'FontName',parm.fontnm,'LineWidth',1.5)
    % set(gcf,'PaperUnits','Inches')
    % set(gcf,'PaperSize',[12/2.57 12/2.57])
    % set(gcf,'color','w')
   
    % pbaspect([1 1 0.5]) 
    % axis tight
     % xmin = 0;
     % xmax = 24;
% %     ymin = -300.0;
% %     ymax = +50.0;
     % xlim([xmin xmax]);
     % xt = [xmin:(xmax-xmin)/12:xmax];
     % xticks(xt);
% %     ylim([ymin ymax]);
% %     yt = [ymin:(ymax-ymin)/7:ymax];
% %     yticks(yt);
% %     %view(46,20);
    % parm.fn = parm.fns(i).name;
    % print(gcf,[parm.fn '.PasModalVsPas.png'],'-dpng',parm.plotresu); 
    % close(gcf);
    

    % %% plot 2d histogram
    % figure('Position',[100 100 1920 1920]);
    % box('on')
    % %% variant 2 - build in histogram2d plot matlab
    % Xedges = [-Inf 0:0.1:24 Inf]; %distance in multiples of cells
    % Yedges = [-Inf 0:0.5:63 Inf]; %disori
    % h = histogram2(GrainIDsDisoriToModal(2,:).*256,GrainIDsDisoriToModal(4,:), ...
        % Xedges,Yedges, ...
        % 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    % colorbar;
    % if contains(parm.fn,'.20')
        % caxis([0.0 0.10])
    % end
    % if contains(parm.fn,'.100')
        % caxis([0.0 0.020])
    % end
    % if contains(parm.fn,'.345')
        % caxis([0.0 0.002])
    % end
    % view(0,90);
    % xlabel({'Distance to boundary'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    % ylabel({'Disorientation to grain modal ori (°)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
    % set(gca,'FontSize',parm.fontsz,'FontName',parm.fontnm,'LineWidth',1.5)
    % set(gcf,'PaperUnits','Inches')
    % set(gcf,'PaperSize',[12/2.57 12/2.57])
    % set(gcf,'color','w')
    % pbaspect([1 1 0.5]) 
    % axis tight
    % xmin = 0;
    % xmax = 24;
    % xlim([xmin xmax]);
    % xt = [xmin:(xmax-xmin)/12:xmax];
    % xticks(xt);
    % parm.fn = parm.fns(i).name;
    % print(gcf,[parm.fn '.ActModalVsAct.png'],'-dpng',parm.plotresu); 
    % close(gcf);

