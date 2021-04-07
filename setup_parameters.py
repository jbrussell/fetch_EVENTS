########## 1-get_event_waveforms ##########
minmagnitude = 7.0 #6.5; 7.0; # Lower magnitude cutoff for event search
webservice = "IRIS"
network = "YO" # YO ENAM

input_stalist = 1 # 0 if use all stations
if input_stalist: # List of stations
    stalist = '/Users/russell/Lamont/ENAM/DATA/stalist_good.txt'
    text_file = open(stalist, "r")
    stations = text_file.read().split('\n')
    text_file.close()
    stations = ','.join(stations).replace(" ", "")
else: # Use all available stations
    stations = "*"
    
tstart = "2014-04-13T00:00:00"
tend = "2015-03-29T00:00:00"
is_downsamp = 1 # Downsample?
sr_new = 1 # Downsample Hz (samples/sec)
trlen = 6000 # Length of traces (sec)
comps = []
comps = ["LHZ", "LH1", "LH2"] #["BHZ", "BH1", "BH2", "BDH"] # Components to download
is_removeresp = 1 # Remove response?
outunits = 'DISP' # DISP, VEL, ACC [For pressure channels, should use "VEL"]
isCMT_params = 1 # Use GCMT parameters for SAC header or default IRIS? (time, lat, lon, depth)
iscentroid = 1 # if isCMT_params=1, use centroid instead of epicentral?
search_dir = './' + webservice + '_' + network + '_' + str(minmagnitude) + '/' # OUTPUT PATH

########## 2-rotate_H1H2_NE_RT ##########
# Read orientations
ori_path = '/Users/russell/Lamont/ENAM/DATA/orientations/YO_orientations.txt'
iscleandir = 0 # Clean directories by removing all BHT, BHR, BHE, and BHN before calculating.
# Define Naming Conventions
Zcomp = "LHZ" # input (same as above)
H1comp = "LH1" # input (same as above)
H2comp = "LH2" # input (same as above)
Ncomp = "LHN" # output
Ecomp = "LHE" # output
Tcomp = "LHT" # output
Rcomp = "LHR" # output


########## 3-CMT2idagrn ##########
CMT2idagrn_path = search_dir+'CMT2idagrn/' # OUTPUT PATH