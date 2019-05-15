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
    parm.fns = dir(['*828*']);
    parm.N = length(parm.fns);
    parm.NX = int32(256);
    parm.NXYZ = int32(parm.NX^3);
        
%% run post-processing loop
for rk=1:1:parm.N
    parm.fn = parm.fns(rk).name;
    display(parm.fn);
    parm.NR = str2double(extractBetween(parm.fn,'NR.','.NC'));
    parm.NC = str2double(extractBetween(parm.fn,'NC.','.bin'));
    parm.ID = str2double(extractBetween(parm.fn,'SimID.','.Rank'));
    parm.RK = str2double(extractBetween(parm.fn,'Rank.','.My'));
        
%% load content
    filename = parm.fn; %'/mnt/home/m.kuehbach/DAMASKPDT_VoroComposer_FINAL/Benchmarking/DAMASKPDT.SimID.30101345.Rank.0.MyProfiling.csv';
    delimiter = ';';
    startRow = 4;
    formatSpec = '%s%f%f%C%C%f%f%f%f%f%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    What = dataArray{:, 1};
    Increment = dataArray{:, 2};
    ID = dataArray{:, 3};
    Category = dataArray{:, 4};
    ParallelismInfo = dataArray{:, 5};
    ProcessVirtualMemory = dataArray{:, 6};
    ProcessResidentSetSize = dataArray{:, 7};
    WallClock = dataArray{:, 8};
    CumulatedWallClock = dataArray{:, 9};
    WallClockFraction = dataArray{:, 10};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;
       
%% load raw heavy data from damaskpdt output
    N = length(What(:,1));
    What = cellstr(What(:,1));
    Category = cellstr(Category(:,1));
    cand = contains(Category(:,1),'IO');
    IORem = sum(WallClock(cand(1:N),1)); %get all I/O stuff
    NonIORem = sum(WallClock(~cand(1:N),1)); %get all non I/O stuff
    Total = sum(WallClock(:,1)); % get total
    
    %get most costly non-io from per rank
    j = N;
    for i=N:-1:1
        if cand(i,1)~=1
            A1C = WallClock(i,1);
            A1W = What(i,1);
            j = i;
            break;
        end
    end
    for i=j-1:-1:1
        if cand(i,1)~=1
            A2C = WallClock(i,1);
            A2W = What(i,1);
            j = i;
            break;
        end
    end
    for i=j-1:-1:1
        if cand(i,1)~=1
            A3C = WallClock(i,1);
            A3W = What(i,1);
            break;
        end
    end
          
   % three most costly non io
   WHAT(1+parm.RK,1:7) = [parm.RK,A1W,A2W,A3W,'NonIORem','IORem','Total'];
   MPI(1+parm.RK,1:7) = [parm.RK,A1C,A2C,A3C,Total-A1C-A2C-A3C-IORem,IORem,Total];
   
   %% main memory consumption at peak
   MEM(1+parm.RK,1:3) = [parm.ID,max(ProcessVirtualMemory),max(ProcessResidentSetSize)];
       
    
   clearvars -except rk parm WHAT MPI MEM;
    
    %fid = fopen(parm.csv, 'wt');
    %fprintf(fid, '%s\n','DistanceBinEndPixel;Counts;0.0;0.01;0.05;0.25;0.5;0.75;0.95;0.99;1.0;KSDecision0.05;KSSign0.05;KSDecision0.01;KSSign0.01');
    %fclose(fid);
    %dlmwrite(parm.csv,QUANTMATRIX,'-append','delimiter',';');
end
plot(MPI(:,1),MPI(:,2),'.');
%method spec wwall clock s
min(MPI(:,2))
max(MPI(:,2))
%virt GB
min(MEM(:,2))
max(MEM(:,2))
%resi GB
min(MEM(:,3))
max(MEM(:,3))




%ID, virtual mem,  resident set size
MEM(:,2:3) = MEM(:,2:3)/(1024^3); %GB

%MEM(:,2:3) = MEM(:,2:3)/ ((256^3)*(8*(1+1+4+9+9+9))); 
%considering additional material points etc. Fe was not loaded %(1024^3);
