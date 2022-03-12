% Run idagrn6_sac_excite to calculate full synthetic seismograms and get
% mode excitation for each branch
%
% Must first run run_mineos_check and mk_kernels
%
% JBR 07/18
clear; close all;

COMP = 'Z'; %'Z' 'R' 'T' % Component
is_overwrite = 0; % overwrite previous calculations?

parameter_FRECHET;
TYPE = param.TYPE;
CARDID = param.CARDID;
% EVTPATH = param.EVTPATH;
% STAPATH = param.STAPATH;
% SYNTH_OUT = param.SYNTH_OUT;

STAPATH = param.STAPATH; % input station path
CMTPATH = param.CMTPATH; % input CMT path
% DATAPATH = '~/BROWN/RESEARCH/PROJ_NoMelt/DATA/EVENTS/fetch_EVENTS/IRIS_ZA_5.5_Zcorr/';

if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
    if COMP == 'Z' || COMP == 'R'
        error('Toroidal mode: Component must be T');
    end
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
    if COMP == 'T'
        error('Spheroidal mode: Component must be Z or R');
    end
end

setpath_plotwk;
setpath_idagrn;

evtpaths = dir([CMTPATH,'evt_*']);
for iev = 1:length(evtpaths)
    EVTPATH = [CMTPATH,evtpaths(iev).name];
    tkn = strsplit(evtpaths(iev).name,'_');
    evid = tkn{2};
    
    SYNTH_OUT = [param.IDAGRN,'SYNTH/',param.CARDID,'_b',num2str(N_modes),'/',evid,'/'];
    
    %% Change environment variables to deal with gfortran
    setenv('GFORTRAN_STDIN_UNIT', '5') 
    setenv('GFORTRAN_STDOUT_UNIT', '6') 
    setenv('GFORTRAN_STDERR_UNIT', '0')

%     SAC_OUT = [SYNTH_OUT,'full'];
    SAC_OUT = [SYNTH_OUT];
    if ~exist(SAC_OUT)
        mkdir(SAC_OUT);
    end
    
    filcheck = dir([SAC_OUT,'*',COMP,'.sac']);
    if ~isempty(filcheck) && ~is_overwrite
        disp(['Already processed ',evid,' skipping...'])
        continue
    end

%     EXCITE_OUT = [SYNTH_OUT,'excitation'];
%     if ~exist(EXCITE_OUT)
%         mkdir(EXCITE_OUT);
%     end

    %% Run plot_wk
    write_plotwk(TYPE,CARDID)

    com = ['cat run_plotwk.',lower(TYPE),' | plot_wk > plot_wk.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at plot_wk')
    end
    
    %% Build station list for input
    STAPATH_MINEOS = './stations_mineos.txt';
    stas = readtable(STAPATH);
    stas.Var4 = stas.Var4/1000;
    Nstas = length(stas.Var4);
    fid = fopen(STAPATH_MINEOS,'w');
    fprintf(fid,'%d,\n',Nstas);
    for ista = 1:Nstas
        fprintf(fid,'%s %12.6f %12.6f %10.4f\n',stas.Var1{ista},stas.Var2(ista),stas.Var3(ista),stas.Var4(ista));
    end
    fclose(fid);
    
    %% Run idagrn6_sac_excite
    write_idagrn(TYPE,CARDID,EVTPATH,STAPATH_MINEOS,LENGTH_HR,DT,COMP)

    fprintf('------- Calculating Full synthetics & Excitation %s0-%d: %s-------\n',TYPE,N_modes-1,COMP)
    system(['cat run_idagrn.',lower(TYPE),' > idagrn.in']);
    com = ['cat idagrn.in | idagrn6_sac_excite > idagrn.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at idagrn6_sac_excite')
    end

    system(sprintf('mv *.%s.sac %s',COMP,SAC_OUT));
    if TYPE == 'T'
%         system(sprintf('mv *.excite.asc %s',EXCITE_OUT));
        system('rm *.excite.asc');
    end
    
%     % Rename files to match fetch_events and add magnitude to structure
%     files = dir([SAC_OUT,'/*.',COMP,'.sac']);
%     for ista = 1:length(files)
%         tkn = strsplit(files(ista).name,'.');
%         STA = tkn{1};
%     end

    %% Change the environment variables back to the way they were
    setenv('GFORTRAN_STDIN_UNIT', '-1') 
    setenv('GFORTRAN_STDOUT_UNIT', '-1') 
    setenv('GFORTRAN_STDERR_UNIT', '-1')

    delete('idagrn.in','idagrn.LOG','plot_wk.LOG',['*.',lower(TYPE)])
end