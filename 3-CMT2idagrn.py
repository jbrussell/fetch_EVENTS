# %% markdown
# # CMT2idagrn
# 
# Reads GCMT events and builds evt_######## files which are input to idagrn to make synthetic seismograms
# 
# ##### JBR - 2/7/18
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
from obspy import read
import numpy as np
import matplotlib.pylab as plt
%matplotlib inline
from datetime import datetime
import os
import calendar
import urllib

# %% codecell
# Load Client
client = Client(webservice)
print(client)

if not os.path.exists(CMT2idagrn_path):
    os.makedirs(CMT2idagrn_path)
# %% codecell
t1 = UTCDateTime(tstart)
t2 = UTCDateTime(tend)

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


# %% codecell
for iev, ev in enumerate(cat_evts) :
    # determine index of desired origin type
    ior = 0 
    if isCMT_params:
        for origin in cat_evts[0].origins:
            if origin.origin_type == ortype:
                break
            ior = ior + 1
    
    # extract date for naming folders
    time = ev.origins[ior].time 
    
    date = datetime.strptime(str(time),'%Y-%m-%dT%H:%M:%S.%fZ')
    evname = date.strftime('%Y%m%d%H%M')
    
    const = 1e19
    scale = 0.1
    lat = ev.origins[ior].latitude
    lon = ev.origins[ior].longitude
    depth = ev.origins[ior].depth
    m_rr = ev.focal_mechanisms[0].moment_tensor.tensor.m_rr/const
    m_tt = ev.focal_mechanisms[0].moment_tensor.tensor.m_tt/const
    m_pp = ev.focal_mechanisms[0].moment_tensor.tensor.m_pp/const
    m_rt = ev.focal_mechanisms[0].moment_tensor.tensor.m_rt/const
    m_rp = ev.focal_mechanisms[0].moment_tensor.tensor.m_rp/const
    m_tp = ev.focal_mechanisms[0].moment_tensor.tensor.m_tp/const
    
    print('Working on '+CMT2idagrn_path+'evt_'+evname)
    f = open(CMT2idagrn_path+'evt_'+evname,'w')
    f.write("%s\n%.3f %.3f %.2f\n%.1f\n%.5f\n%.5f\n%.5f\n%.5f\n%.5f\n%.5f\n" % (evname,lat,lon,depth/1000,scale,m_rr,m_tt,m_pp,m_rt,m_rp,m_tp))
    f.close()
    
    
    
# %% codecell
