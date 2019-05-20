%% DAMASKPDT get number of grains
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
    parm.fns = dir(['*428*GrainData.csv']);
    parm.N = length(parm.fns);
    parm.NX = int32(256);
    parm.NXYZ = int32(parm.NX^3);
        
%% run post-processing loop
for rk=1:1:parm.N
    parm.fn = parm.fns(rk).name;
    display(parm.fn);
    parm.ID = str2double(extractBetween(parm.fn,'SimID.','.Rank'));
    parm.RK = str2double(extractBetween(parm.fn,'Rank.','.Incr'));
    parm.Incr = str2double(extractBetween(parm.fn,'345.Incr.','.GrainData.csv'));
        
%% load content
    filename = parm.fn; %'/mnt/home/m.kuehbach/DAMASKPDT_VoroComposer_FINAL/Benchmarking/DAMASKPDT.SimID.30101345.Rank.0.MyProfiling.csv';
    delimiter = ';';
    startRow = 4;
    formatSpec = '%f%f%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    GrainID = dataArray{:, 1};
    NumberOfIPsAssigned = dataArray{:, 2};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
       
%% load raw heavy data from damaskpdt output
   N = length(unique(GrainID(:,1)));
   %GID = cellstr(GrainID(:,1));
          
   % three most costly non io
   %WHAT(1+parm.RK,1:7) = [parm.RK,A1W,A2W,A3W,'NonIORem','IORem','Total'];
   NGR(1+parm.RK,1) = parm.RK;
   NGR(1+parm.RK,2) = parm.Incr;
   NGR(1+parm.RK,3) = N;
      
   clearvars -except rk parm NGR;
end

fid = fopen(['DAMASKPDT.SimID.' num2str(parm.ID) '.NumberOfGrainsReconstructed.csv'], 'wt');
fprintf(fid, '%s;%s;%s\n','Rank;Increment;NumberOfGrains');
for i=1:1:length(NGR(:,1))
    fprintf(fid, '%d;%d;%d\n',NGR(i,1),NGR(i,2),NGR(i,3));
end
fclose(fid);
%dlmwrite(parm.csv,QUANTMATRIX,'-append','delimiter',';');
%plot(0:1:27,NGR(:,1),'.');