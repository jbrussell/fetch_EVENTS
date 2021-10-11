# %% markdown
# # writeCMTSOLUTION
# 
# Write out standard CMTSOLUTION files for each event
# 
# ##### JBR - 10/11/2021
# %% codecell
# %load_ext autoreload
# %autoreload
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
# %matplotlib inline
from datetime import datetime
import os
import calendar
import urllib

# %% codecell
# Load Client
client = Client(webservice)
print(client)

if not os.path.exists(CMTSOLUTION_path):
    os.makedirs(CMTSOLUTION_path)
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


# %% codecell
for iev, ev in enumerate(cat_evts) :
    
    evname = ev.origins[0].time.strftime('%Y%m%d%H%M');
    outpath = CMTSOLUTION_path+'/'+evname+'_CMTSOLUTION'
    print('Working on '+outpath)
    
    ev.write(outpath,format='CMTSOLUTION')
    
    
