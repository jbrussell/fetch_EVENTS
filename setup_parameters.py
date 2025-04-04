########## 1-get_event_waveforms ##########
# Network parameters
tstart = "2014-04-13T00:00:00"
tend = "2015-03-29T00:00:00"
webservice = "IRIS"
networks = ["YO"] # YO ENAM
minlatitude = 31 # min latitude of stations
maxlatitude = 38 # max latitude of stations
minlongitude = -78 # min longitude of stations
maxlongitude = -70 # max longitude of stations
comps = ["BHZ", "BH1", "BH2", "BDH"] # Components to download

input_stalist = 0 # 0 if use all stations
if input_stalist: # List of stations
    stalist = '/Users/russell/Lamont/ENAM/DATA/stalist_good.txt'
    text_file = open(stalist, "r")
    stations = text_file.read().split('\n')
    text_file.close()
    stations = ','.join(stations).replace(" ", "")
else: # Use all available stations
    stations = "*"

# Event parameters
minmagnitude = 7.0 #6.5; 7.0; # Lower magnitude cutoff for event search
maxmagnitude = 10.0 #6.5; 7.0; # upper magnitude cutoff for event search
mindepth = 0 # [km] minimum event depth
maxdepth = 50 # [km] maximum event depth
isCMT_params = 1 # Use GCMT parameters for SAC header or default IRIS? (time, lat, lon, depth)
iscentroid = 1 # if isCMT_params=1, use centroid instead of epicentral?

# Waveform processing parameters    
is_downsamp = 1 # Downsample?
sr_new = 1 # Downsample Hz (samples/sec)
trlen = 6000 # Length of traces (sec)
# WARNING! List the full channel names. Do not use wildcards. Bad things will happen...
is_removeresp = 1 # Remove response?
outunits = 'DISP' # DISP, VEL, ACC [For pressure channels, should use "VEL"]

search_dir = './' + webservice + '_M' + str(minmagnitude) + '_' + str(mindepth) + '_' + str(maxdepth) + 'km' + '/'

########## 2-rotate_H1H2_NE_RT ##########
# Read H1 orientations relative to North; clockwise positive (only used in 2-rotate-H1H2_NE_RT.py, not needed for 2-rotate_NE_RT.py)
ori_path = '/Users/russell/Lamont/ENAM/DATA/orientations/YO_orientations.txt'
iscleandir = 0 # Clean directories by removing all BHT, BHR, BHE, and BHN before calculating.
# Define Naming Conventions
Zcomp = "BHZ" # input (same as above)
H1comp = "BH1" # input (same as above)
H2comp = "BH2" # input (same as above)
Ncomp = "BHN" # output
Ecomp = "BHE" # output
Tcomp = "BHT" # output
Rcomp = "BHR" # output


# ########## 3-CMT2idagrn ##########
CMT2idagrn_path = search_dir+'CMT2idagrn/' # OUTPUT PATH

# ########## 4-writeCMTSOLUTION ##########
CMTSOLUTION_path = search_dir+'CMTSOLUTIONS/' # OUTPUT PATH
