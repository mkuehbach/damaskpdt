%% DAMASKPDT get benchmarking results
% Markus K\"uhbach, 2019/05/09 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;
    parm.plotresu = '-r600';
    parm.fontsz = 24;
    parm.fontnm = 'DejaVu Sans';
    parm.fontszax = 22;
    
    
%% load user data
    parm.fns = dir(['*GrainStressSpatDistr*.bin']);
    parm.N = length(parm.fns);
    parm.NX = 256;
    parm.NXYZ = parm.NX^3;
    
    
%% run post-processing loop
for i=1:1:parm.N
    parm.fn = parm.fns(i).name;
    display(parm.fn);
    parm.NR = str2double(extractBetween(parm.fn,'NR.','.NC'));
    parm.NC = str2double(extractBetween(parm.fn,'NC.','.bin'));
    parm.ID = str2double(extractBetween(parm.fn,'SimID.','.Rank'));
       
%% load raw heavy data from damaskpdt output
    display(['Loading ' parm.fn ' ...']);
    fid = fopen(parm.fn);
    StressSpatDistr = fread(fid,[parm.NC, parm.NR],'double'); %implicit transpose
	fclose(fid);
    display([parm.fn ' loaded']);
    clearvars -except i cs ss parm StressSpatDistr;
	
	%we dont need to correlate grain specific stresses for now but we could, so for now kill first row with the grain IDs
	StressSpatDistr(1,:) = [];	
	%now distance values cached in StressSpatDistr
	Distances(1,:) = StressSpatDistr(1,:); %previously 2nd row now the first!
    StressSpatDistr(1,:) = [];
    
    % go to MPa
    StressSpatDistr = StressSpatDistr ./ 1.0e6;
    
%% correlate stress tensor values with distance
	% successively walk over 3x3 tensor component values a11, a12, a13, a21, ...
	for c=1:1:9
		
		%% plot 2d histogram
%		figure('Position',[100 100 1920 1920]);
%		box('on')
		
% 		%% plotting core
% 		%if opt.binning == 1
% 		%    %% variant 1 - 2d kernel density
% 		%    [bandwidth,density,X,Y]=kde2d(DAT);
% 		%    surf(X,Y,density,'LineStyle','none'), view([0,90])
% 		%    colormap parula, hold on, alpha(1.0)
% 		%    colorbar
% 		%end
% 		%if opt.binning == 2
% 
% 		%% variant 2 - build in histogram2d plot matlab
% 		Xedges = [-Inf 0:0.2:24 Inf]; %distance in multiples of cells
% 		Yedges = [-Inf -300.0:10.0:+300.0 Inf]; %stress in Pa
% 		h = histogram2(Distances(1,:).*256,StressSpatDistr(c,:), ...
% 			Xedges,Yedges, ...
% 			'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
% 		colorbar;
%         if contains(parm.fn,'.20')
%             caxis([0.0 0.100])
%         end
% 		if contains(parm.fn,'.100')
%             caxis([0.0 0.01])
%         end
%         if contains(parm.fn,'.345')
% 			caxis([0.0 0.010])
% 		end
% 		view(0,90);
% 		xlabel({'Distance to boundary'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 		%if opt.dislotype == 1
% 		switch c
% 			case 1
% 				ylabel({'\sigma_{11} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 2
% 				ylabel({'\sigma_{12} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 3
% 				ylabel({'\sigma_{13} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 4
% 				ylabel({'\sigma_{21} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 5
% 				ylabel({'\sigma_{22} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 6
% 				ylabel({'\sigma_{23} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 7
% 				ylabel({'\sigma_{31} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 8
% 				ylabel({'\sigma_{32} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			case 9
% 				ylabel({'\sigma_{33} (MPa)'},'FontSize',parm.fontsz,'FontName',parm.fontnm);
% 			otherwise
% 		end				
% 		%if opt.dislotype == 2
% 		%ylabel({'log_{10}(\Sigma\rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
% 		%end
% 		set(gca,'FontSize',parm.fontsz,'FontName',parm.fontnm,'LineWidth',1.5)
% 		set(gcf,'PaperUnits','Inches')
% 		set(gcf,'PaperSize',[12/2.57 12/2.57])
% 		set(gcf,'color','w')
% 		pbaspect([1 1 0.5])
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
 		parm.fn = parm.fns(i).name;
% 		switch c
% 			case 1
% 				print(gcf,[parm.fn '.Sigma11.png'],'-dpng',parm.plotresu);
% 			case 2
% 				print(gcf,[parm.fn '.Sigma12.png'],'-dpng',parm.plotresu);
% 			case 3
% 				print(gcf,[parm.fn '.Sigma13.png'],'-dpng',parm.plotresu);
% 			case 4
% 				print(gcf,[parm.fn '.Sigma21.png'],'-dpng',parm.plotresu);
% 			case 5
% 				print(gcf,[parm.fn '.Sigma22.png'],'-dpng',parm.plotresu);
% 			case 6
% 				print(gcf,[parm.fn '.Sigma23.png'],'-dpng',parm.plotresu);
% 			case 7
% 				print(gcf,[parm.fn '.Sigma31.png'],'-dpng',parm.plotresu);
% 			case 8
% 				print(gcf,[parm.fn '.Sigma32.png'],'-dpng',parm.plotresu);
% 			case 9
% 				print(gcf,[parm.fn '.Sigma33.png'],'-dpng',parm.plotresu);
% 		end				
% 		close(gcf);
% 		display([parm.fn ' processed completely']);

        %% compute quantile distributions
        nquants = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95,0.99];
        nq = length(nquants);
        dbins = [0.0:0.2:24.0]./parm.NX; %end values
        Nbins = length(dbins)-1;
        
        %% build reference for KS Test
        dright = 0.0/parm.NX;
        dleft =  1.0/parm.NX;
        %refcand = [];
        %refvalues = [];
        %refcand = Distances(1,:) > dright & Distances(1,:) <= dleft; 
        %refvalues = StressSpatDistr(c,refcand);
   
        QUANTMATRIX = nan(Nbins,1+1+1+nq+1); %+1+1+1+1); %binid, cnts, 
        %min, quant values, max, accept yes/no 0.05, p@0.05, accept yes/no @0.01, p@0.01
        %figure
        for d=1:1:Nbins
            clearvars cand values;
            dright = dbins(d);
            dleft = dbins(d+1);
            QUANTMATRIX(d,1) = dleft * parm.NX;
            cand = Distances(1,:) > dright & Distances(1,:) <= dleft; 
            QUANTMATRIX(d,2) = sum(cand);
            values = StressSpatDistr(c,cand);
            if ~isempty(values)
                QUANTMATRIX(d,2+1) = min(values);
                QUANTMATRIX(d,2+1+1:2+1+nq) = quantile(values,nquants);
                QUANTMATRIX(d,2+1+nq+1) = max(values);
                %end
                %hold on
                %ecdf(values);
                % KS test significance for difference of the distributions
                %clearvars h5 p5 h1 p1;
                %[h5, p5] = kstest2(values,refvalues,'Alpha',0.05); %[h,p] if h == true rejecting null hypothesis at the 0.01 significance level
                %[h1, p1] = kstest2(values,refvalues,'Alpha',0.01);
                %QUANTMATRIX(d,2+1+nq+1+1) = double(h5); %acce y/n
                %QUANTMATRIX(d,2+1+nq+1+2) =  p5; %significance @0.05
                %QUANTMATRIX(d,2+1+nq+1+3) = double(h1); %acce y/n
                %QUANTMATRIX(d,2+1+nq+1+4) =  p1; %significance @0.05
            end
            display(d);
        end
        display(['Statistical analysis for component ' num2str(c) ' done']);
%         figure
%         hold on
%         plot(QUANTMATRIX(:,1),QUANTMATRIX(:,2+1),'.');
%         for q=1:1:nq
%             hold on
%             plot(QUANTMATRIX(:,1),QUANTMATRIX(:,2+1+q),'.');
%         end
%         hold on
%         plot(QUANTMATRIX(:,1),QUANTMATRIX(:,2+1+nq+1),'.');
%         figure
%         semilogy(QUANTMATRIX(:,1),QUANTMATRIX(:,2),'.');
        switch c
			case 1
				parm.csv = [parm.fn '.Sigma11.csv'];
			case 2
				parm.csv = [parm.fn '.Sigma12.csv'];
			case 3
				parm.csv = [parm.fn '.Sigma13.csv'];
			case 4
				parm.csv = [parm.fn '.Sigma21.csv'];
			case 5
				parm.csv = [parm.fn '.Sigma22.csv'];
			case 6
				parm.csv = [parm.fn '.Sigma23.csv'];
			case 7
				parm.csv = [parm.fn '.Sigma31.csv'];
			case 8
				parm.csv = [parm.fn '.Sigma32.csv'];
			case 9
				parm.csv = [parm.fn '.Sigma33.csv'];
        end
        fid = fopen(parm.csv, 'wt');
        fprintf(fid, '%s\n','DistanceBinEndPixel;Counts;Min;0.01;0.05;0.25;0.5;0.75;0.95;0.99;Max'); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
        for d=1:1:Nbins
            fprintf(fid, '%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f;%.18f\n',QUANTMATRIX(d,:)); %KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
        end
        fclose(fid);
        %dlmwrite(parm.csv,QUANTMATRIX,'-append','delimiter',';');
    end    
	clearvars -except i parm; 
end