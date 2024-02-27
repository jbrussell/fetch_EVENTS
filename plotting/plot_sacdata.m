%% PLOT SEISMOGRAMS

clear; close all;
%plot native
fig1 = figure(1); clf;
set(gcf,'position',[141    28   938   633]);
%figure(2); clf;

data_dir = 'IRIS_YQ_6.5'; %'IRIS_XX_5.5_Zcorr'; %'IRIS_ZA_7.0';
comp = 'BHZ';
fb_min = 1/150; % 150 sec
fb_max = 1/20; % 20 sec

% data_dir = 'IRIS_YQ_6.5'; %'IRIS_XX_5.5_Zcorr'; %'IRIS_ZA_7.0';
% comp = 'BHT';
% fb_min = 1/150; %1/100; % 150 sec
% fb_max = 1/20; % 20 sec

trace_space = 0; %0.15; % spacing between traces
amp = 0.3 / 2;
max_dist_deg = 150;

groupv = 4.3;
winrefmin = -1200; %-800;
winrefmax = 1500; %1000;

figpath = ['./figs/',data_dir,'/',num2str(1/fb_max),'_',num2str(1/fb_min),'/'];
if ~exist(figpath)
    mkdir(figpath);
end

%% LOAD EVENT LIST
evs = dir(['./',data_dir,'/*']);
% evs_pat = dir('/Users/russell/Lamont/nomelt_DATA/GSDF_LOVE_cleanEVTs/sacdata_T/2012*');
num_evs = size(evs,1);

%% LOAD SAC FILES
land = shaperead('landareas', 'UseGeoCoords', true);
load coast;
for j = 1:num_evs
    figure(1); clf
    sac_fils = dir(['./',data_dir,'/',evs(j).name,'/*',comp,'.sac']);
%     sac_fils_pat = dir(['/Users/russell/Lamont/nomelt_DATA/GSDF_LOVE_cleanEVTs/sacdata_T/',evs(j).name,'/*.sac']);
    if isempty(sac_fils)
        continue
    end
    
    num_fil = size(sac_fils,1);
    DIST_min = 9999;
    DIST_max = 0;
    SAC = {};
    dists = 0;
    data = [];
    itrace = 0;

    for i = 1:num_fil
        PATH_sac_fils = ['./',data_dir,'/',evs(j).name,'/',sac_fils(i).name];
%         display("Checking:")
%         display(sac_filsT(i).name);
%         display(sac_filsZ(i).name);
        SAC{i} = rdsac(PATH_sac_fils);
        if i == 1
            tref = round(SAC{i}.HEADER.DIST/1000/groupv);
            twindow = [tref+winrefmin:tref+winrefmax];
            twindow = twindow(twindow>=0);
        end
        
        if SAC{i}.HEADER.GCARC > DIST_max
            DIST_max = SAC{i}.HEADER.GCARC;
        end
        if SAC{i}.HEADER.GCARC < DIST_min
            DIST_min = SAC{i}.HEADER.GCARC;
        end
        EVDP = SAC{i}.HEADER.EVDP/1000;
        MAG = SAC{i}.HEADER.MAG;
        
        %Filter data
        d_xt{i} = SAC{i}.d;
        %Taper data
        d_xt{i} = d_xt{i}.*tukeywin(length(d_xt{i}),0.08);
%         ts = SAC{i}.t;
        ts = 0:length(d_xt{i})-1;
        fs = 1/(ts(2)-ts(1));
        [b,a] = butter(2,[fb_min/(fs/2) fb_max/(fs/2)]); % (20 - 150 seconds)
        %fvtool(b,a);
        d_xt{i} = filtfilt(b,a,d_xt{i});
        
        GCARC = SAC{i}.HEADER.GCARC;
        
        figure(1); box off;
        %Plot data
        if GCARC > max_dist_deg
            continue
        end
        if min(abs(GCARC-dists)) > trace_space
            itrace = itrace+1;
            dists(itrace,:) = GCARC;
            dwindow = d_xt{i}(twindow+1);
            data(itrace,:) = amp*(dwindow/max(abs(dwindow)));
        end
%         box on;
            
    end
    if isempty(data)
        continue
    end
    figure(1); box off;
    plot(twindow,data+dists,'linewidth',1.5,'color',[0 0 0]); hold on;
    xlabel('Time (s)');
    ylabel('Distance (\circ)');
    title(['Event:',evs(j).name,'  ',num2str(EVDP),' km  M',num2str(MAG)]);
    set(gca,'fontsize',16,'linewidth',1);
%         xlim([ts(1) ts(floor(end/1.9))]);
    xlim([twindow(1) twindow(end)]);
    ylim([DIST_min-0.2 DIST_max+0.25]);
    
    % Plot inset map
    axes('Position',[.75 .79 .35*.6 .4*.6])
    box off;
%         axesm('robinson','MapLatLimit',[-90 90],'maplonlimit',[45,44])
    h = axesm('robinson','MapLatLimit',[-90 90],'maplonlimit',[-180,+180]+SAC{1}.HEADER.STLO);
    setm(gca,'grid','off','frame', 'on');
    axis off
    p = findobj(h,'type','patch'); % Find background
    set(p,'FaceColor',[1 1 1]); % Change background to white
%         geoshow(land,"FaceColor",[0.7 0.7 0.7])
    plotm(lat, long,'k');
%         raypath_lat = [SAC{1}.HEADER.STLA,SAC{1}.HEADER.EVLA];
%         raypath_lon = [SAC{1}.HEADER.STLO,SAC{1}.HEADER.EVLO];
    [raypath_lat,raypath_lon] = gcwaypts(SAC{1}.HEADER.STLA,SAC{1}.HEADER.STLO,SAC{1}.HEADER.EVLA,SAC{1}.HEADER.EVLO,50);
    plotm(raypath_lat,raypath_lon,'-r','linewidth',1.5);
    plotm(SAC{1}.HEADER.STLA,SAC{1}.HEADER.STLO,'ok','linewidth',0.5,'markerfacecolor','b','markersize',10);
    plotm(SAC{1}.HEADER.EVLA,SAC{1}.HEADER.EVLO,'pk','linewidth',0.5,'markerfacecolor','r','markersize',20);
        

    save2pdf([figpath,'/',comp,'_evt_',evs(j).name,'_M',num2str(MAG),'_',num2str(1/fb_max),'_',num2str(1/fb_min),'s','.pdf'],fig1,1000);
end