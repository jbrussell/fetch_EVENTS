# %% markdown
# # Download Event Data
# Searches for events greater than given magnitude and downloads specified components, resamples to specified rate, removes instrument response, and saves as SAC file.
# 
# ##### JBR - 2/3/18
# ##### JBR - 2/7/18 : Include option to use GCMT parameters in SAC headers
# ##### JBR - 3/10/21: Use IRIS SPUD to get GCMT parameters rather than the GCMT webpage
# %% codecell
%load_ext autoreload
%autoreload
from setup_parameters import *
import matplotlib.pyplot as plt
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy.io.sac import SACTrace
from obspy.geodetics import gps2dist_azimuth, locations2degrees
import numpy as np
import os
from datetime import datetime
import calendar
import urllib

# %% codecell
if not os.path.exists(search_dir):
    os.makedirs(search_dir)
    
# LOAD CLIENT
client = Client(webservice)
print(client)
# %% codecell
# LOAD EVENT CATALOGUE
t1 = UTCDateTime(tstart)
t2 = UTCDateTime(tend)
if isCMT_params==1 :
    # Load events from GCMT catalogue using IRIS SPUD
    url_query = 'http://ds.iris.edu/spudservice/momenttensor/ids?' \
               +'evtstartdate='+t1.strftime('%Y-%m-%dT%H:%M:%S') \
               +'&evtenddate='+t2.strftime('%Y-%m-%dT%H:%M:%S') \
               +'&evtminmag='+str(minmagnitude)
    evids = urllib.request.urlopen(url_query)
    events_str = '&'.join([line.decode("utf-8").replace("\n", "") for line in evids])+'&'
    url_ndk = 'http://ds.iris.edu/spudservice/momenttensor/bundleids/ndk?'+events_str
    cat_evts = obspy.read_events(url_ndk)
    if iscentroid: # Use centroid parameters
        ortype = 'centroid'
    else: # Use hypocenter parameters
        ortype = 'hypocenter'
else: # Read events from IRIS catalogue
    cat_evts = client.get_events(starttime=t1, endtime=t2, minmagnitude=minmagnitude)

# %% codecell
# LOAD STATIONS
inventory = client.get_stations(network=network, station=stations, channel=comps[0], starttime=t1, endtime=t2)
inventory.plot(projection="local",label=False)
fig = inventory.plot(show=False) 
cat_evts.plot(fig=fig)  
fig.savefig(search_dir+"events.pdf", bbox_inches="tight")

file = open('stations.txt', 'w')
for ista in range(0,len(inventory[0])) :
    file.write("%5s %12f %12f %12f\n" % (inventory[0].stations[ista]._code, 
                                        inventory[0].stations[ista]._latitude, 
                                        inventory[0].stations[ista]._longitude, 
                                        inventory[0].stations[ista]._elevation))
file.close()


# %% codecell
# Get event times
f = open(search_dir+'evlist.txt','w')
for iev, ev in enumerate(cat_evts) :  
    
    # determine index of desired origin type
    ior = 0 
    if isCMT_params:
        for origin in cat_evts[0].origins:
            if origin.origin_type == ortype:
                break
            ior = ior + 1
    
    # Get earthquake parameters
    tbeg = ev.origins[ior].time
    tend = tbeg + trlen
    evdp = ev.origins[ior].depth
    evla = ev.origins[ior].latitude
    evlo = ev.origins[ior].longitude 
    mag = ev.magnitudes[0].mag
    
    # extract date for naming folders
    date = datetime.strptime(str(tbeg),'%Y-%m-%dT%H:%M:%S.%fZ')
    evname = date.strftime('%Y%m%d%H%M')
    evdir = search_dir + evname + '/'
    f.write(evname+'\n')
    if not os.path.exists(evdir):
        os.makedirs(evdir)
        print('Working on : ',evname)
    else:
        print('Skipping '+evname+'...')
        continue
        
    # Loop through stations
    for ista in range(0,len(inventory[0])) :
        stel = inventory[0].stations[ista].elevation
        stla = inventory[0].stations[ista].latitude
        stlo = inventory[0].stations[ista].longitude
        station = inventory[0].stations[ista].code
        vals = gps2dist_azimuth(lat1=stla, lon1=stlo, lat2=evla, lon2=evlo)
        dist = vals[0]
        baz = vals[1]
        az = vals[2]
        gcarc = locations2degrees(lat1=stla, long1=stlo, lat2=evla, long2=evlo)
    
        # Loop through components
        for icomp, comp in enumerate(comps):   
            try:
                st = client.get_waveforms(network=network, station=station, location="*", channel=comp, starttime=tbeg, endtime=tend, attach_response=True)
            except:
                print('Missing data for station: ',station)
                continue
            if len(st) > 1: # Check for data gaps and fill with 0's
                st.merge(method=1, fill_value=0)
            sr = st[0].stats.sampling_rate
            if is_removeresp:
                try:
                    # Check whether pressure channel, if so use "VEL" option which doesn't add or remove zeros
                    if st[0].stats.response.instrument_sensitivity.input_units == 'PA':
                        st.remove_response(output='VEL', zero_mean=True, taper=True, taper_fraction=0.05, pre_filt=[0.001, 0.005, sr/3, sr/2], water_level=600)
                    else:
                        st.remove_response(output=outunits, zero_mean=True, taper=True, taper_fraction=0.05, pre_filt=[0.001, 0.005, sr/3, sr/2], water_level=600)
                except:
                    print('Skipping... issue removing response: '+station+' '+evname)
                    continue
            st.trim(starttime=tbeg, endtime=tend, pad=True, nearest_sample=False, fill_value=0) # make sure correct length
            st.detrend(type='demean')
            st.detrend(type='linear')
            st.taper(type="cosine",max_percentage=0.05)
#             sr_old = st[0].stats.sampling_rate
            if is_downsamp==1:
                st.filter('lowpass', freq=0.4*sr_new, zerophase=True) # anti-alias filter
                st.resample(sampling_rate=sr_new)
                st.detrend(type='demean')
                st.detrend(type='linear')
                st.taper(type="cosine",max_percentage=0.05)
    
            # convert to SAC and fill out station/event header info
            sac = SACTrace.from_obspy_trace(st[0])
            sac.stel = stel
            sac.stla = stla
            sac.stlo = stlo
            sac.evdp = evdp
            sac.evla = evla
            sac.evlo = evlo
            sac.mag = mag
            sac.dist = dist
            sac.az = az
            sac.baz = baz
            sac.gcarc = gcarc
            
            sac_out = evdir + evname + '.' + network + '.' + station + '.' + comp + '.sac'
            sac.write(sac_out)
f.close()
# %% codecell

# %% codecell
