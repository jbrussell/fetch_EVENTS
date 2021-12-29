# %% markdown
# # Rotate To R-T
# Reads in Z-H1-H2 and first rotates to Z-N-E, then to Z-R-T and write out sac file. A left hand convention is assumed where H2 is 90$^{\circ}$ clockwise from H1, and Z points up. Input instrument azimuths are H1 clockwise of north.
# 
# ### Station orientation file is of the form:
# STA1 H1-AZI1
# 
# STA2 H1-AZI2
# 
# STA3 H1-AZI3
# 
# ...
# 
# If station is missing from list of orientations, then it is simply skipped.
# 
# ##### JBR - 2/3/18
# %% codecell
%load_ext autoreload
%autoreload
from setup_parameters import *
import numpy as np
from obspy import read
from obspy.signal.rotate import rotate2zne, rotate_ne_rt
import pandas as pd
import matplotlib.pylab as plt
from obspy.io.sac import SACTrace
import os
%matplotlib inline
# %% codecell
data_dir = search_dir

# Read events
evs = pd.read_csv(data_dir+'/evlist.txt',' ',header=None,skipinitialspace=True)
evs.columns = ["event"]

# Orientations
ori = pd.read_csv(ori_path,' ',header=None,skipinitialspace=True)
ori.columns = ["sta", "H1azi"]
stas = ori['sta']
ori = ori.set_index('sta')
# %% codecell
for iev, event in enumerate(evs['event']):
    event = str(event)
    evdir = data_dir + event + '/'
    print('Working on : ',event)
    if iscleandir:
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Tcomp+'.sac')
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Rcomp+'.sac')
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Ncomp+'.sac')
        os.system('rm '+evdir + event+'.'+'*'+'.'+'*'+'.'+Ecomp+'.sac')
    
    for ista, sta in enumerate(stas):
        try:
            st = read(evdir + event+'.*.'+str(sta)+'.*.sac', debug_headers=True)
        except Exception:
            print('Missing data for station: ',str(sta))
            continue
        H1azi = ori.loc[sta]['H1azi']
        for itr in range(0,len(st)):
            if st[itr].stats.channel == H1comp:
                h1 = st[itr].data
            elif st[itr].stats.channel == H2comp:
                h2 = st[itr].data
            elif st[itr].stats.channel == Zcomp:
                z = st[itr].data
        ba = st[0].stats.sac.baz
        
        # Rotate Z-H1-H2 to Z-N-E
        traces_zne = rotate2zne(data_1=z , azimuth_1=0, dip_1=-90, 
                   data_2=h1, azimuth_2=H1azi, dip_2=0, 
                   data_3=h2, azimuth_3=H1azi+90, dip_3=0)
        z2 = traces_zne[0]
        n = traces_zne[1]
        e = traces_zne[2]
        
        # Rotate N-E to R-T
        traces_rt = rotate_ne_rt(n=n, e=e, ba=ba)
        r = traces_rt[0]
        t = traces_rt[1]

        # Define new data streams
        st_bhn = st[0].copy()
        st_bhn.stats.channel = Ncomp
        st_bhn.data = n
        st_bhe = st[0].copy()
        st_bhe.stats.channel = Ecomp
        st_bhe.data = e
        st_bht = st[0].copy()
        st_bht.stats.channel = Tcomp
        st_bht.data = t
        st_bhr = st[0].copy()
        st_bhr.stats.channel = Rcomp
        st_bhr.data = r

        # Remove existing file
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Ncomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Ncomp+'.sac')
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Ecomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Ecomp+'.sac')
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Tcomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Tcomp+'.sac')
        if os.path.exists(evdir + event+'.'+network+'.'+str(sta)+'.'+Rcomp+'.sac'):
            os.remove(evdir + event+'.'+network+'.'+str(sta)+'.'+Rcomp+'.sac')
        
        network = st_bhr.stats.network
        # Save BHN, BHE, BHR, and BHT
        sac_n = SACTrace.from_obspy_trace(st_bhn)
        sac_n.write(evdir + event+'.'+network+'.'+str(sta)+'.'+Ncomp+'.sac')
        sac_e = SACTrace.from_obspy_trace(st_bhe)
        sac_e.write(evdir + event+'.'+network+'.'+str(sta)+'.'+Ecomp+'.sac')
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
