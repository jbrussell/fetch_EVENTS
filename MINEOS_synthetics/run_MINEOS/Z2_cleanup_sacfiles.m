% Copy header structure from sac data files to make synthetics look
% identical to the data
%
% JBR 06/2021
clear; close all;

parameter_FRECHET;
TYPE = param.TYPE;
CARDID = param.CARDID;
COMP = param.COMP; %'Z' 'R' 'T' % Component
COMP_prefix = param.COMP_prefix; % channel prefix
LENGTH_HR = param.LENGTH_HR;
DT = param.DT;
% EVTPATH = param.EVTPATH;
% STAPATH = param.STAPATH;
% SYNTH_OUT = param.SYNTH_OUT;

% Data directory (to copy header info)
DATAPATH = param.DATAPATH;
% New output path for synthetic data
NEW_SYNTH_OUT = param.NEW_SYNTH_OUT;


SYNTH_PATHS = [param.IDAGRN,'SYNTH/',param.CARDID,'_b',num2str(N_modes),'/'];
evtpaths = dir(SYNTH_PATHS); evtpaths=evtpaths(~ismember({evtpaths.name},{'.','..'}));
for iev = 1:length(evtpaths)
    evid = evtpaths(iev).name;    
    sacfiles = dir([SYNTH_PATHS,evid,'/*.',COMP,'.sac']);
    if isempty(sacfiles)
        continue
    end
    
    for ista = 1:length(sacfiles)
        tkn = strsplit(sacfiles(ista).name,'.');
        STA = tkn{1};
        
        % Load synthetic
        fullpathsac = [SYNTH_PATHS,evid,'/',sacfiles(ista).name];
        if ~exist(fullpathsac)
            continue
        end
        synth = rdsac(fullpathsac);
        
        % Load data
        datafile = dir([DATAPATH,'/',evid,'/',evid,'*.',STA,'.*',COMP_prefix,COMP,'.sac']);
        if isempty(datafile)
            disp(['Skip missing sac file in data ',fullpathsac]);
            continue
        end
        if length(datafile)>1
            error('More than one matching data file... check component');
        end
        data = rdsac([DATAPATH,'/',evid,'/',datafile.name]);
        
        % Replace syntheric header values with data values
        flds = fields(data.HEADER);
        for ifld = 1:length(flds)
            fld = flds{ifld};
            if strcmp(fld,'DELTA') || strcmp(fld,'NPTS')
                continue
            end
            synth.HEADER.(fld) = data.HEADER.(fld);
        end
        synth.HEADER.KSTNM = strtrim(synth.HEADER.KSTNM);
        synth.HEADER.KCMPNM = strtrim(synth.HEADER.KCMPNM);
        synth.HEADER.KNETWK = strtrim(synth.HEADER.KNETWK);
        
        % Write out new sac file
        H = synth.HEADER;
        sacdata = synth.d;
        opath = [NEW_SYNTH_OUT,'/',evid,'/'];
        if ~exist(opath)
            mkdir(opath)
        end
        startDate = datetime(H.NZYEAR,1,H.NZJDAY); % use this to convert jday to month - day
        startYear = num2str(year(startDate), '%04d');
        startMonth = num2str(month(startDate), '%02d');
        startDay = num2str(day(startDate), '%02d');
        fullevid = [startYear,startMonth,startDay,num2str(H.NZHOUR,'%02d'),num2str(H.NZMIN,'%02d'),num2str(H.NZSEC,'%02d'),num2str(H.NZMSEC,'%03d')];
        startTime = datenum(fullevid,'yyyymmddhhMMSSFFF');
        sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.sac',opath,evid, H.KNETWK, H.KSTNM, H.KCMPNM));
        mksac(sac_path,sacdata,startTime,H);        

    end

end