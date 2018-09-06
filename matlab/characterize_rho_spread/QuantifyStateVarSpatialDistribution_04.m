%% damaskpdt
% functionality to quantify the spatial distribution of state variable
% values with respect to the proximity at local interfaces
% Markus Kühbach, 2018/05/26 m.kuehbach at mpie.de
    clear;
    clc;
    digits 32;
    format long;

%% user interaction
    dbl.simid = 10236;
    dbl.INCR = [0,10,20,30,40,50,60,70,80,90,100,114,128,142,156,170,184,198,212,226,240,265,290,315,340,365,390,415];
    dbl.prefix = ['E:/Paper13_OrientationGradients/cc_analysis/1000/' num2str(dbl.simid) '/'];
    dbl.NXYZ = 128^3;
    dbl.h = 1.0/((dbl.NXYZ)^(1/3));
    dbl.NC = 13;

%% analyses options
    opt.binning = 2;    %1 - ks2d, 2 - hist2d (binning no ksdensity!)
    opt.dislotype = 1;  %1 - edges, 2 - dipoles, either case \sum all slipsystems per data point
    if opt.dislotype == 1
        opt.dtype = 'EdgeDnsSpatDistr';
        dbl.fns = dir([dbl.prefix '*.' opt.dtype '*bin']);
    else opt.dislotype == 2
        opt.dtype = 'DipoleDnsSpatDistr';
        dbl.fns = dir([dbl.prefix '*.' opt.dtype '*bin']);
    end
    
%% batch processing
for i=1:1:1 %2:1:length(dbl.INCR(1,:))
   clearvars -except dbl opt i SIGNIFICANCE;
%% load heavy data
    %there is no natural sort on dbl.fns(:).name therefore find corresponding index for INCR(i)
    for k=1:1:length(dbl.fns)
        if opt.dislotype == 1
            thisone = str2double(extractBetween(dbl.fns(k).name,'0.1.415.Incr.', '.Edge'));
        end
        if opt.dislotype == 2
            thisone = str2double(extractBetween(dbl.fns(k).name,'0.1.415.Incr.', '.Dipole'));
        end
                
        if thisone == dbl.INCR(i)
            thisone = k;
            break;
        end
    end
    parms.fn_f64 = [dbl.prefix dbl.fns(k).name];
    NC = dbl.NC;
    NR = str2double(extractBetween(parms.fn_f64,'NR.','.NC'));
    fid = fopen(parms.fn_f64);
    DAT = fread(fid,[NC,NR],'double'); %implicit transpose
    fclose(fid);
    display(['Loaded ' parms.fn_f64]);


%% debugging heavy data
    %     clear;
    %     clc;
    %     digits 32;
    %     format long;
    %     %fid = fopen(['X:\damaskpdt\run\GrainReconDebug\' ...
    %     %    'DAMASKPDT.SimID.1001.Incr.0.EdgeDnsSpatDistr.NR.4728.NC.13.bin']);
    %     fid = fopen(['X:\damaskpdt\run\GrainReconDebug\' ...
    %         'DAMASKPDT.SimID.1001.Incr.0.DipoleDnsSpatDistr.NR.4728.NC.13.bin']);
    % 
    %     DAT1 = fread(fid,[13,4728],'double'); %implicit transpose
    %     fclose(fid);
    %     %fid = fopen(['X:\damaskpdt\run\GrainReconDebug\' ...
    %     %    'DAMASKPDT.SimID.1002.Incr.0.EdgeDnsSpatDistr.NR.4728.NC.13.bin']);
    %     fid = fopen(['X:\damaskpdt\run\GrainReconDebug\' ...
    %         'DAMASKPDT.SimID.1002.Incr.0.DipoleDnsSpatDistr.NR.4728.NC.13.bin']);
    %     DAT2 = fread(fid,[13,4728],'double'); %implicit transpose
    %     fclose(fid);
    %     for i=1:13
    %         clearvars TMP;
    %         for j=1:4728
    %             TMP(1,j) = abs(DAT1(i,j)-DAT2(i,j));
    %         end
    %         sum(TMP(1,:))
    %         clearvars TMP;
    %     end
    

%% compute rho sum
    RHO(1,:) = DAT(1,:); %distances to interface
    RHO(2,:) = DAT(2,:);
    for c=3:NC
        RHO(2,:) = RHO(2,:) + DAT(c,:);
    end

%% descriptive stuff
    %min(RHO(1,:))
    %max(RHO(1,:))


%% plot the distribution
    figure('Position',[100 100 2000 1000]);
    %grid('on')
    box('on')
    %view(0,90)
    %view(58,24)
    fontszcap = 22;
    fontszax = 22;
    fontnm = 'Calibri';

%% plotting core
    if opt.binning == 1
        %% variant 1 - 2d kernel density
        [bandwidth,density,X,Y]=kde2d(DAT);
        surf(X,Y,density,'LineStyle','none'), view([0,90])
        colormap parula, hold on, alpha(1.0)
        colorbar
    end
    if opt.binning == 2
        %% for dislocation densities rather plot log10 values as otherwise
        % bin length in each dimension differ by order of magnitudes

        %% variant 2 - build in histogram2d plot matlab
        %h = histogram2(DAT(:,1),DAT(:,2),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');
        %h = histogram2(RHO(1,:)./dbl.h,log10(RHO(2,:)),'XBinEdges',[1:1:15].*dbl.h, 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
        h = histogram2(RHO(1,:)./dbl.h,log10(RHO(2,:)),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
        colorbar;
    end

    xlabel({'Distance to boundary (\cdot h)'},'FontSize',fontszcap,'FontName',fontnm);
    if opt.dislotype == 1
        ylabel({'log_{10}(\Sigma\rho^s_{edge}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
    end
    if opt.dislotype == 2
        ylabel({'log_{10}(\Sigma\rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
    end

    set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
    set(gcf,'PaperUnits','Inches')
    set(gcf,'PaperSize',[30 30])
    set(gcf,'color','w')
    pbaspect([1 2 0.5])
    xmin = 0.0;
    xmax = 14.0;
    % ymin = min(RHO(:,2));
    % ymax = max(RHO(:,2));
    xlim([xmin xmax]);
    xt = [xmin:(xmax-xmin)/14:xmax];
    xticks(xt);
    % ylim([ymin ymax]);
    % yt = [ymin:(ymax-ymin)/5:ymax];
    % yticks(yt);
    view(90,-90);
    %view(46,20);
    
    %very evil code, xlabels and ylabels does times disappear even though
    %visible, Matlab plotting is simply fucked up when compared to
    %AfterEffects etc...
    %hax = get(gca,'XLabel');
    %set(hax,'Visible','on')
    %set(hax,'Position',[5.65 12.417 -0.062],'visible','on');
    %hay = get(gca,'YLabel');
    %set(hay,'Position',[16.585 13.07 -0.099]),'visible','on');
    if opt.dislotype == 1
        print(gcf,[parms.fn_f64 '.png'],'-dpng','-r500'); close(gcf);
    end
    if opt.dislotype == 2
        print(gcf,[parms.fn_f64 '.png'],'-dpng','-r500'); close(gcf);
    end
%clearvars -except dbl opt;

    
%% compute normalized mono-parameter eCDF for \rho for binning in distance classes with spacing equivalent to the discretization
    NX = dbl.NXYZ^(1/3);
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
        display([num2str(b) '--> [' num2str(binleft) '  ' num2str(binright) ' ] ' num2str(length(cand(1,:))/dbl.NXYZ)]);
        %display(median(log10(cand)))
    
        if ( ~isempty(cand) )
            %h --> 0 do not reject, h --> 1 reject null hypothesis
            [h, p] = kstest2(log10(cand),ref(1,:),'Alpha',alpha); %[h,p] if h == true rejecting null hypothesis at the 0.01 significance level
            SIGNIFICANCE(b,1) = dbl.INCR(i);
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
%             SIGNIFICANCE(b,1) = dbl.INCR(i);
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