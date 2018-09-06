%% damaskpdt
% functionality to quantify the spatial distribution of state variable
% values with respect to the proximity at local interfaces
% Markus Kühbach, 2018/02/02 m.kuehbach at mpie.de
clear;
clc;
digits 32;
format long;

%% user interaction
%%20180426 --- fn = 'Z:\damaskpdt\run\500g256x256x256\DAMASKPDT.SimID.250.Incr.250.DipoleDensitySpatDistr.csv';
%%20180426 --- DAT = load_statevarspatdistr2(fn);
%fn = 'X:\damaskpdt\run\500g256x256x256\DAMASKPDT.SimID.320.Incr.320.EdgeDensitySpatDistr.csv';
%fn = 'Z:/damaskpdt/run/500g128x128x128_compZ/DAMASKPDT.SimID.1002.Incr.415.EdgeDensitySpatDistr.csv';

%analyses options
opt_binning = 2; %1 - ks2d, 2 - hist2d (binning no ksdensity!)
opt_dislotype = 1; %1 - edges, 2 - dipoles, either case \sum all slipsystems per data point

%input file details
SimID = 1002;
IncrID = 415;
NR = 128^3;
NC = 13;
if opt_dislotype == 1
    dtype = 'EdgeDnsSpatDistr';
else opt_dislotype == 2
    dtype = 'DipoleDnsSpatDistr';
end
mydir = 'Z:/damaskpdt/run/500g128x128x128_compZ/';
prefix  = ['DAMASKPDT.SimID.' num2str(SimID) '.Incr.' num2str(IncrID) '.' dtype];
suffix = ['NR.' num2str(NR) '.NC.' num2str(NC) '.bin'];
%fn = ['X:\damaskpdt\run\GrainReconDebug\DAMASKPDT.SimID.0.Incr.0.EdgeDnsSpatDistr' suffix];
fn = [mydir prefix '.' suffix];


%% load heavy data
display(['Loading ' fn ' ...']);
fid = fopen(fn);
DAT = fread(fid,[NC,NR],'double'); %implicit transpose
fclose(fid);


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
for c=3:13 %check against NC
    RHO(2,:) = RHO(2,:) + DAT(c,:);
end

%% descriptive stuff start with the total dislocation density
min(RHO(1,:))
max(RHO(1,:))

%% compute normalized mono-parameter eCDF for \rho for binning in distance classes with spacing equivalent to the discretization
NX = 128;
SubStep = 2.0;
clearvars i binleft binright cand h p SIGNIFICANCE;
inc = 2;
for i=1:12
    binleft = (i-1)/(SubStep*NX); %assuming implicit scaling
    binright = i/(SubStep*NX);
    cand = RHO(2, RHO(1,:) >= binleft & RHO(1,:) < binright);
    if i == 1
        candref = log10(cand);
    end
    
    %display([num2str(i) '-->' num2str(binleft) '  ' num2str(binright) '  ' num2str(length(cand(1,:)))]);
    display(median(cand))
    
    if ( length(cand(1,:)) > 0 )
        [h, p] = kstest2(log10(cand),candref,'Alpha',0.05); %[h,p] if h == true rejecting null hypothesis at the 0.01 significance level
        SIGNIFICANCE(inc).Reject(i) = h;
        SIGNIFICANCE(inc).PValue(i) = p;
        SIGNIFICANCE(inc).Median(i) = median(cand);
    else
        SIGNIFICANCE(inc).Reject(i) = true;
        SIGNIFICANCE(inc).PValue(i) = NaN;
        SIGNIFICANCE(inc).Median(i) = NaN;
    end
    clearvars cand;
end
%two-sided KolmogorovSmirnov reject against null hypothesis data are drawn
%from the same, non-parametric assumption free

%% plot the evolution of median and hypothesis testing results against strain
targetclass = 1;
% ##start building the figure
figure('Position',[100 100 1000 1000]);
%grid('on')
box('on')
%view(0,90)
fontszcap = 22;
fontszax = 22;
fontnm = 'Calibri';
for inc=1:length(SIGNIFICANCE(:))
    % different classes different color
    hold on
    plot(0.1,log10(SIGNIFICANCE(inc).Median(targetclass)),'.');
    %plot(SIGNIFICANCE(inc).Strain(targetclass),SIGNIFICANCE(inc).Median(targetclass),'.'); 
end
set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
set(gcf,'PaperUnits','Inches')
set(gcf,'PaperSize',[30 30])
set(gcf,'color','w')
pbaspect([2 2 2])
xmin = 0.0;
xmax = 0.32;
ymin = 15;
ymax = 16;
xlim([xmin xmax]);
xt = [xmin:(xmax-xmin)/8:xmax];
xticks(xt);
ylim([ymin ymax]);
yt = [ymin:(ymax-ymin)/10:ymax];
yticks(yt);

xlabel({'\epsilon_{vM}'},'FontSize',fontszcap,'FontName',fontnm);
if opt_dislotype == 1
    ylabel({'log_{10}(\Sigma_{s} \rho^s_{edge}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
end
if opt_dislotype == 2
    ylabel({'log_{10}(\Sigma_{s} \rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
end

if opt_dislotype == 1
    print(gcf,[fn '.EdgeTotal.png'],'-dpng','-r500'); close(gcf);
end
if opt_dislotype == 2
    print(gcf,[fn '.DipoleTotal.png'],'-dpng','-r500'); close(gcf);
end


    




%% plot the distribution
figure('Position',[100 100 2000 1000]);
%grid('on')
box('on')
view(0,90)
view(58,24)
fontszcap = 22;
fontszax = 22;
fontnm = 'Calibri';

if opt_binning == 1
    %% variant 1 - 2d kernel density
    [bandwidth,density,X,Y]=kde2d(DAT);
    surf(X,Y,density,'LineStyle','none'), view([0,90])
    colormap parula, hold on, alpha(1.0)
    colorbar
elseif opt_binning == 2
    %% for dislocation densities rather plot log10 values as otherwise
    % bin length in each dimension differ by order of magnitudes

    %% variant 2 - build in histogram2d plot matlab
    %h = histogram2(DAT(:,1),DAT(:,2),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');
    h = histogram2(RHO(1,:),log10(RHO(2,:)),'Normalization','probability','FaceColor','flat'); %'ShowEmptyBins','on');  
    colorbar;
else
    disp('Choose a plotting mode!');
end

xlabel({'Normalized distance'},'FontSize',fontszcap,'FontName',fontnm);
if opt_dislotype == 1
    ylabel({'log_{10}(\Sigma_{s} \rho^s_{edge}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
end
if opt_dislotype == 2
    ylabel({'log_{10}(\Sigma_{s} \rho^s_{dipol}) (m^{-2})'},'FontSize',fontszcap,'FontName',fontnm);
end

set(gca,'FontSize',fontszcap,'FontName',fontnm,'LineWidth',1.5)
set(gcf,'PaperUnits','Inches')
set(gcf,'PaperSize',[30 30])
set(gcf,'color','w')
pbaspect([1 2 0.5])
% xmin = min(DAT(:,1));
% xmax = max(DAT(:,1));
% ymin = min(RHO(:,2));
% ymax = max(RHO(:,2));
% xlim([xmin xmax]);
% xt = [xmin:(xmax-xmin)/5:xmax];
% xticks(xt);
% ylim([ymin ymax]);
% yt = [ymin:(ymax-ymin)/5:ymax];
% yticks(yt);
if opt_dislotype == 1
    print(gcf,[fn '.EdgeTotal.png'],'-dpng','-r500'); close(gcf);
end
if opt_dislotype == 2
    print(gcf,[fn '.DipoleTotal.png'],'-dpng','-r500'); close(gcf);
end


















zlim([0 0.005])


modelFun =  @(p,x) 1.0/(2*pi.*(p(1)^2)) .* exp( -((x-p(2))./(2.0*p(1))).^2 );
startingVals = [0.5 14.5];
coefEsts = nlinfit(h.YBinEdges(1:225),h.Values(1,:), modelFun, startingVals);
xgrid = linspace(0,20,100);
line(xgrid, modelFun(coefEsts, xgrid), 'Color','r');


% % clear all
% % N = 100;
% % X=rand(N,1);
% % Y=rand(N,1); %sin(X*10*pi)+randn(size(X))/3; 
% % data=[X,Y];
% % 
% % plot(data(:,1),data(:,2),'.')
% % min(data)
% % max(data)
% % 
% % % apply routine
% % [bandwidth,density,X,Y]=kde2d(data);
% % % plot the data and the density estimate
% % surf(X,Y,density,'LineStyle','none'), view([0,90])
% % colormap parula, hold on, alpha(1.0)
% % colorbar
% % 
% % % show data points
% % set(gca, 'color', 'blue');
% % plot(data(:,1),data(:,2),'w.','MarkerSize',5)