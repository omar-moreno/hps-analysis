#
#
#
#
#
#
#

#--- Imports ---#
#---------------#
import argparse
import os
import sys

from ROOT import TFile, TCanvas, TH1F, AddressOf, TH2F, TGraphErrors
from ROOT import TMultiGraph, gSystem 
from ROOT import RooRealVar, RooDataHist, RooPlot, RooArgSet, RooArgList
from ROOT import RooGaussian, RooFit, RooAddPdf, RooFFTConvPdf, RooGaussModel, RooAddModel

#--- Main ---#
#------------#

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="DST file to process.")
args = parser.parse_args()

# If a DST file was not specified, exit
if not args.input: 
	print "\nPlease specify a DST file to process."
	print "\nUse the -h flag for usage.\n"
	sys.exit(2)

# Open the DST file
dst_file = TFile(str(args.input))
if not dst_file.IsOpen():
	print "\nDST file " + str(args.input) + " couldn't be opened!"
	sys.exit(2)

# Get the TTree
dst_tree = dst_file.Get("HPS_Event")

# Load the DST libraries
hps_dst_path = os.environ['HPS_DST_HOME']
hps_dst_path += "/build/lib/libHpsEvent.so"
gSystem.Load(hps_dst_path)

from ROOT import HpsEvent

hps_event = HpsEvent()
b_hps_event = dst_tree.GetBranch("Event")
b_hps_event.SetAddress(AddressOf(hps_event))

canvas = TCanvas("canvas", "Extrapolation Analysis", 800, 800)

h_track_d0 = TH1F("h_track_d0", "Track D0", 100, -2, 2)
h_track_d0.SetLineColor(3)
h_gbl_track_d0 = TH1F("h_gbl_track_d0", "GBL Track D0", 100, -2, 2)
h_gbl_track_d0.SetLineColor(2)

h_sp_residuals_bp_vs_p_5hit = TH2F("h_sp_residuals_bp_vs_p_5hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_sp_residuals_bp_vs_p_6hit = TH2F("h_sp_residuals_bp_vs_p_6hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_sp_residuals_nbp_vs_p_5hit = TH2F("h_sp_residuals_nbp_vs_p_5hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_sp_residuals_nbp_vs_p_6hit = TH2F("h_sp_residuals_nbp_vs_p_6hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 160, -20, 20)

track_parameters = [0]*5
gbl_track_parameters = [0]*5

for entry in xrange(0, dst_tree.GetEntries()):
	
	# Print the event number every 500 events
	if (entry+1)%500 == 0: print "Processing event " + str(entry+1)
	
	# Get the event from the tree
	dst_tree.GetEntry(entry)

	# Loop over all the Ecal scoring plane hits in the event and match
	# them to their respective GBL track using the original (seed) track
	for ecal_sp_hit_n in xrange(0, hps_event.getNumberOfEcalScoringPlaneHits()):

		ecal_sp_hit = hps_event.getEcalScoringPlaneHit(ecal_sp_hit_n)

		for gbl_track_n in xrange(0, hps_event.getNumberOfGblTracks()):

			gbl_track = hps_event.getGblTrack(gbl_track_n)

			#if(gbl_track.getSeedChi2() != track.getChi2()): continue
			
			#h_track_d0.Fill(gbl_track.getSeedD0())
			#h_gbl_track_d0.Fill(gbl_track.getD0())

			#gbl_track_parameters[0] = gbl_track.getD0()	
			#gbl_track_parameters[1] = gbl_track.getPhi()
			#gbl_track_parameters[2] = gbl_track.getKappa()
			#gbl_track_parameters[3] = gbl_track.getZ0()
			#gbl_track_parameters[4] = math.tan(math.pi/2 - gbl_track.getTheta())

			#track_parameters[0] = track.getD0()
			#track_parameters[1] = track.getPhi()
			#track_parameters[2] = track.getOmega()
			#track_parameters[3] = track.getZ0()
			#track_parameters[4] = track.getTanLambda()

			#extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(gbl_track_parameters)
			#gbl_track_position = extrapolator.extrapolate_track(ecal_sp_hit.getPosition()[2])

			#delta_x = ecal_sp_hit.getPosition()[0] - gbl_track_position[0]
			#delta_y = ecal_sp_hit.getPosition()[1] - gbl_track_position[1]
			
			#extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(track_parameters)
			#track_position = extrapolator.extrapolate_track(ecal_sp_hit.getPosition[2])

			#if sp_n_tracker_hits == 5:
			#	h_sp_residuals_bp_vs_p_5hit.Fill(sp_momentum, delta_x)
			#	h_sp_residuals_nbp_vs_p_5hit.Fill(sp_momentum, delta_y)
			#elif sp_n_tracker_hits == 6:
			#	h_sp_residuals_bp_vs_p_6hit.Fill(sp_momentum, delta_x)
			#	h_sp_residuals_nbp_vs_p_6hit.Fill(sp_momentum, delta_y)


canvas.Print("track_parameters.pdf[")

#--- Fit the D0 plots ---#
#------------------------#

track_d0 = RooRealVar("track_d0", "D0 (mm)", -2, 2)
track_d0_mean = RooRealVar("track_d0_mean", "Track D0 Mean", 0, -2, 2)
track_d0_sigma = RooRealVar("track_d0_sigma", "Track D0 Sigma", .1, 0, 2)
track_d0_gaussian = RooGaussModel("track_d0_gaussian", "Track D0 Gaussian", track_d0, track_d0_mean, track_d0_sigma)

track_d0_data = RooDataHist("track_d0_data", "Track D0 Data", RooArgList(track_d0), h_track_d0)
track_d0_plot = track_d0.frame()
track_d0_data.plotOn(track_d0_plot)

track_d0_gaussian.fitTo(track_d0_data)
track_d0_gaussian.plotOn(track_d0_plot)

track_d0_plot.Draw()
canvas.Print("track_parameters.pdf(")

gbl_track_d0 = RooRealVar("gbl_track_d0", "D0 (mm)", -2, 2)
gbl_track_d0_mean = RooRealVar("gbl_track_d0_mean", "Track D0 Mean", 0, -2, 2)
gbl_track_d0_sigma = RooRealVar("gbl_track_d0_sigma", "Track D0 Sigma", .1, 0, 2)
gbl_track_d0_gaussian = RooGaussModel("gbl_track_d0_gaussian", "Track D0 Gaussian", gbl_track_d0, gbl_track_d0_mean, gbl_track_d0_sigma)

gbl_track_d0_data = RooDataHist("track_d0_data", "Track D0 Data", RooArgList(gbl_track_d0), h_gbl_track_d0)
gbl_track_d0_plot = gbl_track_d0.frame()
gbl_track_d0_data.plotOn(gbl_track_d0_plot)

gbl_track_d0_gaussian.fitTo(gbl_track_d0_data)
gbl_track_d0_gaussian.plotOn(gbl_track_d0_plot)

gbl_track_d0_plot.Draw()
canvas.Print("track_parameters.pdf(")

canvas.Print("track_parameters.pdf]")

