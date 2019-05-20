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
            %if opt.binning == 1
            %    %% variant 1 - 2d kernel density
            %    [bandwidth,density,X,Y]=kde2d(DAT);
            %    surf(X,Y,density,'LineStyle','none'), view([0,90])
            %    colormap parula, hold on, alpha(1.0)
            %    colorbar
            %end
            %if opt.binning == 2
        
            %ecdf(DAT(10,:))
            %% variant 2 - build in histogram2d plot matlab
            %h = histogram2(DAT(:,1),DAT(:,2),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');
            %h = histogram2(RHO(1,:)./parm.h,log10(RHO(2,:)),'XBinEdges',[1:1:15].*parm.h, 'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
            h = histogram2(DAT(1,:),DAT(10,:)./1.0e6,'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
            colorbar;
            caxis([0.0 0.02])
            view(0,90);
            xlabel({'Normal distance to boundary'},'FontSize',fontszcap,'FontName',fontnm);
            %if opt.dislotype == 1
            ylabel({'\sigma_{33} (MPa)'},'FontSize',fontszcap,'FontName',fontnm);
            %end
            %if opt.dislotype == 2
            %ylabel({'log_{10}(\Sigma\rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
            %end
            set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
            set(gcf,'PaperUnits','Inches')
            set(gcf,'PaperSize',[30 30])
            set(gcf,'color','w')
            pbaspect([2 1 0.5])
            xmin = 0.0;
            xmax = 0.12;
            ymin = -300.0;
            ymax = +50.0;
            xlim([xmin xmax]);
            xt = [xmin:(xmax-xmin)/12:xmax];
            xticks(xt);
            ylim([ymin ymax]);
            yt = [ymin:(ymax-ymin)/7:ymax];
            yticks(yt);
            %view(46,20);

            %very evil code, xlabels and ylabels does times disappear even though
            %visible, Matlab plotting is simply fucked up when compared to
            %AfterEffects etc...
            %hax = get(gca,'XLabel');
            %set(hax,'Visible','on')
            %set(hax,'Position',[5.65 12.417 -0.062],'visible','on');
            %hay = get(gca,'YLabel');
            %set(hay,'Position',[16.585 13.07 -0.099]),'visible','on');
            print(gcf,[parm.prefix fn '.png'],'-dpng','-r600'); 
            close(gcf);
            
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