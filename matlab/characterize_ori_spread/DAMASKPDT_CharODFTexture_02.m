%% characterize orientation spread supplemental to damaskpdt via MTex v5.0.3, we utilize SI units
% Markus K\"uhbach, 2018/05/07 m.kuehbach at mpie.de
clear;
clc;
digits(32);
format long;

%% load user data
    %global dbl;
    %global parms;
    
    dbl.prefix = 'Z:/damaskpdt/run/500g128x128x128_compZ/';
    dbl.fns = dir([dbl.prefix '*.QuatCloud1.bin']); %'*.RA.12.RB.16.*.bin']);
    dbl.n = length(dbl.fns);
%% run post-processing loop
    parms.SimID = 1001;
    
for i=23:1:23 %dbl.n:-1:1
    fn = dbl.fns(i).name;
    parms.IncrID = str2double(extractBetween(fn,'Incr.','.Quat'));
       
%% input file details
    %parms.SimID = 0;
    %parms.IncrID = 0;
    parms.NXYZ = 128^3;
    parms.mydir = dbl.prefix; %'X:\damaskpdt\run\500g128x128x128_compZ\';
    parms.prefix  = ['DAMASKPDT.SimID.' num2str(parms.SimID) '.Incr.' num2str(parms.IncrID)];
    parms.fn_u32 = [parms.mydir parms.prefix '.' 'QuatCloud1.bin'];
    parms.fn_f32 = [parms.mydir parms.prefix '.' 'QuatCloud2.bin'];

%% load raw heavy data from damaskpdt output
    display(['Loading ' parms.fn_u32 ' ...']);
    fid = fopen(parms.fn_u32);
    GrainIDs = int32( fread(fid,[1,parms.NXYZ],'uint32') ); %implicit transpose
    fclose(fid);
    
    display(['Loading ' parms.fn_f32 ' ...']);
    fid = fopen(parms.fn_f32);
    Quats = fread(fid,[4,parms.NXYZ],'float32'); %implicit transpose
    fclose(fid);
    
    clearvars -except dbl parms i GrainIDs Quats;

%% customize_mtex
    parms.fontsz = 26;
    parms.fontnm = 'Calibri';
    setMTEXpref('showCoordinates','on');
    setMTEXpref('FontSize',parms.fontsz);
    setMTEXpref('FontName',parms.fontnm);
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
    
%% build orientation list fast
    ori_in = orientation(quaternion(Quats(1:4,:)),cs,ss);
    clearvars Quats;
    
%% compute ODF
    clearvars odf_apprx;
    parms.odf_resu = 10.0*degree;
    odf_apprx = calcODF(ori_in,'resolution',parms.odf_resu);
    
%% characterize ODF
    ['Texture index ' num2str(textureindex(odf_apprx))]
    parms.vol_scatter = 10.0*degree;
    global volfrac;
    volfrac.cube = volume(odf_apprx, orientation('Euler',0.0*degree,0.0*degree,0.0*degree,cs,ss), parms.vol_scatter)*100.0;
    % add further
    % #####
    
%% write report file
    parms.fn_csv = [parms.prefix '.ODFVolFractions.csv'];
    fid = fopen(parms.fn_csv, 'wt');
    fprintf(fid, '%s\n',['Component;Scatter;VolumeFraction']);
    fprintf(fid, '%s\n',[';degree;percent']);
    fprintf(fid, '%s\n',['Component;Scatter;VolumeFraction']);
    fprintf(fid,'%s;%e;%e\n','Cube', parms.vol_scatter/degree, volfrac.cube );
    % report further
    fclose(fid);
    
%% plot ODF classical or reduced, always triclinic
    opt_odfplot = 2; %1 - key section, 2 - classical 18, parula coloring
    if opt_odfplot == 1
        % only key sections
        plot(odf_apprx,'phi2',[0.0 35.0 45.0]*degree,'contourf','silent');
        colorbar
        parms.fn_odf = [parms.mydir parms.prefix '.ODFDefResu.' num2str(parms.odf_resu/degree,'%.0f') '.KeyPhi2Sections.png'];
        parms.fn_odf = [parms.prefix '.ODFDefResu.' num2str(parms.odf_resu/degree,'%.0f') '.KeyPhi2Sections.png'];
    end
    if opt_odfplot == 2
        % classical 5 degree sections
        plot(odf_apprx,'phi2','sections',18,'contourf','silent');
        colorbar
        parms.fn_odf = [parms.prefix '.ODFDefResu.' num2str(parms.odf_resu/degree,'%.0f') '.ClassicalPhi2.png'];
    end
    parms.plotresu = '-r600';
    print(gcf,parms.fn_odf,'-dpng',parms.plotresu); close(gcf);
end