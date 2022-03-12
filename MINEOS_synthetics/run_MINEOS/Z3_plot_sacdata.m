%% PLOT SEISMOGRAMS

clear all; 
fig1 = figure(1); clf;
% set(gcf,'position',[0    0.9167   17.3056    8.0556]);
set(gcf,'position',[141    28   938   633]);
%figure(2); clf;

data_dir = 'IRIS_ZA_5.5_Zcorr_MINEOS/pa5_5km_b2/';
data_dir2 = 'IRIS_ZA_5.5_Zcorr_MINEOS/pa5_5km_b1/';
comp1 = 'BHZ';
comp2 = 'BHZ';
lgd = {'Fund + 1st','Fund only'};


minMag = 5.5;
maxDepth = 50;

% trace_space = km2deg(30); % km/(km/deg) = deg
trace_space = 0; % km/(km/deg) = deg
amp = 0.3;

winrefmin = -1200; %-800;
winrefmax = 1500; %1000;

fb_min = 1/100; %1/150; %1/100; % 150 sec
fb_max = 1/15; %1/60; %1/20; % 20 sec  

% fb_min = 1/100; %1/150; %1/100; % 150 sec
% fb_max = 1/50; %1/60; %1/20; % 20 sec  

figpath = ['./figs/',data_dir,'/',num2str(1/fb_max),'_',num2str(1/fb_min),'/'];
if ~exist(figpath)
    mkdir(figpath);
end

%% LOAD EVENT LIST
evs = dir(['./',data_dir,'/20*']);
% evs_pat = dir('/Users/russell/Lamont/nomelt_DATA/GSDF_LOVE_cleanEVTs/sacdata_T/2012*');
num_evs = size(evs,1);

%% LOAD SAC FILES
for j = 1:num_evs
    sac_filsT = dir(['./',data_dir,'/',evs(j).name,'/*',comp1,'.sac']);
    sac_filsZ = dir(['./',data_dir2,'/',evs(j).name,'/*',comp2,'.sac']);
%     sac_fils_pat = dir(['/Users/russell/Lamont/nomelt_DATA/GSDF_LOVE_cleanEVTs/sacdata_T/',evs(j).name,'/*.sac']);
        
    num_fil = size(sac_filsZ,1);
    
    if num_fil == 0
        continue
    end
    
    DIST_min = 9999;
    DIST_max = 0;
    SAC = {};
    dists = 0;
    data1 = [];
    data2 = [];
    itrace = 0;
    
    ifplotev = 1;
    for i = 1:num_fil
        PATH_sac_filsT = ['./',data_dir,'/',evs(j).name,'/',sac_filsT(i).name];
        PATH_sac_filsZ = ['./',data_dir2,'/',evs(j).name,'/',sac_filsZ(i).name];
%         display("Checking:")
%         display(sac_filsT(i).name);
%         display(sac_filsZ(i).name);
        SAC{i} = rdsac(PATH_sac_filsT);
        SACZ{i} = rdsac(PATH_sac_filsZ);
        
        if SAC{i}.HEADER.MAG < minMag || SAC{i}.HEADER.EVDP/1000 > maxDepth
            ifplotev = 0;
            continue
        end
        if i == 1
            tref = round(SAC{i}.HEADER.DIST/1000/4.3);
            twindow = [tref+winrefmin:tref+winrefmax];
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
        d_xtZ{i} = SACZ{i}.d;
%         ts = SAC{i}.t;
        ts = 0:length(d_xt{i})-1;
        fs = 1/(ts(2)-ts(1));
        [b,a] = butter(2,[fb_min/(fs/2) fb_max/(fs/2)]); % (20 - 150 seconds)
        %fvtool(b,a);
        d_xt{i} = filter(b,a,d_xt{i});
        d_xtZ{i} = filter(b,a,d_xtZ{i});
        
        GCARC = SAC{i}.HEADER.GCARC;
        
        if min(abs(GCARC-dists)) > trace_space
            itrace = itrace+1;
            dists(itrace,:) = GCARC;
            data1(itrace,:) = amp*(d_xt{i}/max(abs(d_xt{i})));
            data2(itrace,:) = amp*(d_xtZ{i}/max(abs(d_xt{i})));
        end
    end
    if ~ifplotev
        continue
    end
    
    figure(1); box off;
    %Plot data
    h2 = plot(ts,data2+dists,'linewidth',1.5,'color',[0 0 0]); hold on;
    h1 = plot(ts,data1+dists,'linewidth',1.5,'color',[1 0 0]); hold on;
    xlabel('Time');
    ylabel('\Delta (deg)');
    title(['Event:',evs(j).name,'  Depth:',num2str(EVDP),' km  M',num2str(MAG),],'fontsize',16);
    legend([h1(1),h2(1)],lgd,'location','eastoutside')
    
    
%         xlim([ts(1) ts(floor(end/1.9))]);
        xlim([twindow(1) twindow(end)]);
        ylim([DIST_min-0.5 DIST_max+0.5]);
        axes('Position',[.75 .1 .35*.6 .4*.6])
        box on;
        axesm('robinson','MapLatLimit',[-90 90],'maplonlimit',[45,44])
        load coast;
        plotm(lat, long,'k');
        plotm(SAC{1}.HEADER.STLA,SAC{1}.HEADER.STLO,'^r','linewidth',5);
        plotm(SAC{1}.HEADER.EVLA,SAC{1}.HEADER.EVLO,'og','linewidth',5);
        raypath_lat = [SAC{1}.HEADER.STLA,SAC{1}.HEADER.EVLA];
        raypath_lon = [SAC{1}.HEADER.STLO,SAC{1}.HEADER.EVLO];
        plotm(raypath_lat,raypath_lon,'-','linewidth',1.5);
        setm(gca,'grid','on','frame', 'on');
        
% %         figure(2); clf;
% %         %worldmap('world');
% %         axesm('robinson','MapLatLimit',[-30 50],'maplonlimit',[-180,-30])
% %         load coast;
% %         plotm(lat, long,'k');
% %         plotm(SAC{1}.HEADER.STLA,SAC{1}.HEADER.STLO,'^r','linewidth',5);
% %         plotm(SAC{1}.HEADER.EVLA,SAC{1}.HEADER.EVLO,'og','linewidth',5);
% %         raypath_lat = [SAC{1}.HEADER.STLA,SAC{1}.HEADER.EVLA];
% %         raypath_lon = [SAC{1}.HEADER.STLO,SAC{1}.HEADER.EVLO];
% %         plotm(raypath_lat,raypath_lon,'-');
% %         setm(gca,'grid','on','meridianlabel','on','frame', 'on','parallellabel','on')
% %     pause;

    save2pdf([figpath,'/',comp1,'_',comp2,'_evt_',evs(j).name,'_M',num2str(MAG),'_',num2str(1/fb_max),'_',num2str(1/fb_min),'s','.pdf'],fig1,1000);
    figure(1); clf
end