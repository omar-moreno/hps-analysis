
import argparse
import sys
import os

import TrackExtrapolator

from PyCintex import gSystem, gDirectory, AddressOf

#--- Main ---#
#------------#

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",  help="ROOT file to process")
parser.add_argument("-o", "--output",  help="Name of output pdf file")
args = parser.parse_args()

# If an input ROOT file was not specified, exit
if not args.input:
	print "\nPlease specify a ROOT file to process."
	print "\nUse the -h flag for usage\n"
	sys.exit(2)

# Load the HpsEvent library
if os.getenv('HPS_DST_HOME') is None: 
	print "Error! Environmental variable HPS_DST_HOME is not set."
	print "\nExiting ..."
	sys.exit(2)
	
hps_dst_path = os.environ['HPS_DST_HOME']
hps_dst_path += "/build/lib/libHpsEvent.so"
gSystem.Load(hps_dst_path)

# import the modules used by HpsEvent i.e. HpsEvent, 
# SvtTrack, EcalCluster ...
from ROOT import HpsEvent, SvtTrack

# Open the ROOT file
root_file = TFile(str(args.input))

# Get the TTree "HPS_EVENT" containing the HpsEvent branch and all
# other colletions
tree = root_file.Get("HPS_Event")

# Create an HpsEvent object in order to read the TClonesArray 
# collections
hps_event = HpsEvent()

# Get the HpsEvent branch from the TTree 
b_hps_event = tree.GetBranch("Event")
b_hps_event.SetAddress(AddressOf(hps_event))

track_parameters = [0]*5

# Loop over all events in the file
for entry in xrange(0, tree.GetEntries()) : 
    
    # Print the event number every 500 events
    if (entry+1)%500 == 0 : print "Event " + str(entry+1)

    # Read the ith entry from the tree.  This "fills" HpsEvent and allows 
    # access to all collections
    tree.GetEntry(entry)
    
    # Loop over all tracks in the event
    for track_n in xrange(0, hps_event.getNumberOfTracks()) : 
        
        # Get the track from the event
        track = hps_event.getTrack(track_n)
       
        # Set the track parameters
        track_parameters[0] = track.getD0()
        track_parameters[1] = track.Phi()
        track_parameters[2] = track.getOmega()
        track_parameters[3] = track.getZ0()
        track_parameters[4] = track.getTanLambda()

        print "Track parameters: " + str(track_parameters)
    
        extrapolator = TrackExtrapolator.TrackExtrapolator(track_parameters)
        
        print "D0: " + str(extrapolator.get_d0())
        print "Phi0: " + str(extrapolator.get_phi0())
        print "x0: " + str(extrapolator.get_x0())

