%% DAMASKPDT pull grain-specific complexity related quantities from the tessellator stdout
% MTex v5.0.3, we utilize SI units
% Markus K\"uhbach, 2019/05/19 m.kuehbach at mpie.de
    clear;
    clc;
    digits(32);
    format long;
    parm.plotresu = '-r600';
    parm.fontsz = 26;
    parm.fontnm = 'Calibri';
    parm.fontszax = 22;

%% load user data
    parm.fn = 'D:\DAMASKPDT_VoroComposer_FINAL\Postprocessing\Benchmarking\MPITotalTime.23489.txt';
    parm.ID = str2double(extractBetween(parm.fn,'MPITotalTime.','.txt'));
    delimiter = {' ',';'};
    formatSpec = '%*s%f%*s%f%*s%*s%f%*s%f%f%f%f%f%f%[^\n\r]';
    fileID = fopen(parm.fn,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false); 
    fclose(fileID);
    TMP = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    
    RESULTS = unique(TMP(:,1));
    NINCR = length(RESULTS(:,1));
    RES(:,1) = reshape(1:1:NINCR,[NINCR,1]);
    RES(:,2) = RESULTS(:,1);
    RES(:,3) = 0; % find all unique grains
    RES(:,4) = 0.0; %total time taken to process all grains in this increment
    RES(:,5) = 0; %total number of facet coverage and edge and vertex tests taken
    for i=1:1:NINCR
        thisone = RES(i,2);
        cand = TMP(:,1) == thisone;
        RES(i,3) = length(unique(TMP(cand,2)));
        RES(i,4) = sum(TMP(cand,3));
        RES(i,5) = sum(TMP(cand,7)) + sum(TMP(cand,8));
    end
    RES(:,6) = RES(:,5) ./ min(RES(:,5));
    clearvars -except parm RES;
    plot(RES(:,1),RES(:,5),'.');

%% write results
    fid = fopen(['DAMASKPDT.SimID.' num2str(parm.ID) '.NumberOfFacetTest.csv'], 'wt');
    fprintf(fid, '%s;%s;%s;%s;%s;%s\n','Incr;DamaskIncrement;NumberOfGrains;TotalTimeExtractContour;SumFacetAndEdgeTests;RelComplexity');
    for i=1:1:length(RES(:,1))
        fprintf(fid, '%.18f;%.18f;%.18f;%.18f;%.18f;%.18f\n', ...
            RES(i,1),RES(i,2),RES(i,3),RES(i,4),RES(i,5),RES(i,6) );
    end
    fclose(fid);