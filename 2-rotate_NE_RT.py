# %% markdown
# # Rotate To R-T
# Reads in Z-N-E and rotates to Z-R-T and write out sac file.
# 
# ##### JBR - 1/30/23
# %% codecell
%load_ext autoreload
%autoreload
from setup_parameters import *
import numpy as np
from obspy import read
from obspy.signal.rotate import rotate_ne_rt
import pandas as pd
import matplotlib.pylab as plt
from obspy.io.sac import SACTrace
import os
%matplotlib inline
# %% codecell
data_dir = search_dir

# Read events
evs = pd.read_csv(data_dir+'/evlist.txt',delimiter=' ',header=None,skipinitialspace=True)
evs.columns = ["event"]

# %% codecell
for iev, event in enumerate(evs['event']):
    event = str(event)
    evdir = data_dir + event + '/'
    print('Working on : ',event)
    if iscleandir:
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Tcomp+'.sac')
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Rcomp+'.sac')
    
    for ista, sta in enumerate(stas):
        try:
            st = read(evdir + event+'.*.'+str(sta)+'.*.sac', debug_headers=True)
        except Exception:
            print('Missing data for station: ',str(sta))
            continue
        for itr in range(0,len(st)):
            if st[itr].stats.channel == Ncomp:
                n = st[itr].data
            elif st[itr].stats.channel == Ecomp:
                e = st[itr].data
            elif st[itr].stats.channel == Zcomp:
                z = st[itr].data
        ba = st[0].stats.sac.baz
        
        # Rotate N-E to R-T
        traces_rt = rotate_ne_rt(n=n, e=e, ba=ba)
        r = traces_rt[0]
        t = traces_rt[1]

        # Define new data streams
        st_bht = st[0].copy()
        st_bht.stats.channel = Tcomp
        st_bht.data = t
        st_bhr = st[0].copy()
        st_bhr.stats.channel = Rcomp
        st_bhr.data = r

        # Remove existing file
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Tcomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Tcomp+'.sac')
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Rcomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Rcomp+'.sac')
        
        network = st_bhr.stats.network
        # Save BHR, BHT
        sac_t = SACTrace.from_obspy_trace(st_bht)
        sac_t.write(evdir + event+'.'+network+'.'+str(sta)+'.'+Tcomp+'.sac')
        sac_r = SACTrace.from_obspy_trace(st_bhr)
        sac_r.write(evdir + event+'.'+network+'.'+str(sta)+'.'+Rcomp+'.sac')
        
#         fmin = 1/100
#         fmax = 1/20
#         st_bht.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
#         st_bhr.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
#         plt.figure(figsize=(10,5))
#         plt.plot(np.arange(0,len(r)), st_bht.data, color="red")
#         plt.plot(np.arange(0,len(t)), st_bhr.data, color="black")
#         # plt.plot(np.arange(0,len(stLHZ[0].data)), stLHZ[0].data, color="blue")
#         plt.xlim(0, 5000)
# #         plt.pause(3)
#         plt.draw()
# %% codecell

# %% codecell
