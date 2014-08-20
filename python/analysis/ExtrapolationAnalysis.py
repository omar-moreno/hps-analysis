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
import math

from ROOT import TFile, TCanvas, TH1F, AddressOf, TH2F, TGraphErrors
from ROOT import TMultiGraph, gSystem 
from ROOT import RooRealVar, RooDataHist, RooPlot, RooArgSet, RooArgList 
from ROOT import RooGaussian, RooFit, RooAddPdf, RooFFTConvPdf, RooGaussModel, RooAddModel
from ROOT import kBlue, kRed, kAzure, kGreen
from ROOT import TLegend

import HpsTrackExtrapolator

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

h_track_d0 = TH1F("h_track_d0", "Track D0", 100, -1.5, 1.5)
h_track_d0.SetLineColor(kAzure+2)
h_track_d0.SetMarkerColor(kAzure+2)
h_gbl_track_d0 = TH1F("h_gbl_track_d0", "GBL Track D0", 100, -1.5, 1.5)
h_gbl_track_d0.SetLineColor(kRed+1)
h_gbl_track_d0.SetMarkerColor(kRed+1)

legend = TLegend(0.7, 0.8, 0.9, 0.9, "", "brNDC")
legend.SetBorderSize(0)
legend.SetFillStyle(0)

legend.AddEntry(h_track_d0, "Seed Tracker", "lp")
legend.AddEntry(h_gbl_track_d0, "GBL Track", "lp")

h_gbl_target_residuals_bp_vs_p_5hit = TH2F("h_gbl_target_residuals_bp_vs_p_5hit", 
		"Bend Plane Residuals at the Target - 5 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_gbl_target_residuals_bp_vs_p_6hit = TH2F("h_gbl_target_residuals_bp_vs_p_6hit", 
		"Bend Plane Residuals at the Target - 6 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_gbl_sp_residuals_bp_vs_p_5hit = TH2F("h_gbl_sp_residuals_bp_vs_p_5hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_gbl_sp_residuals_bp_vs_p_6hit = TH2F("h_gbl_p_residuals_bp_vs_p_6hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_gbl_sp_residuals_nbp_vs_p_5hit = TH2F("h_gbl_sp_residuals_nbp_vs_p_5hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_gbl_sp_residuals_nbp_vs_p_6hit = TH2F("h_gbl_sp_residuals_nbp_vs_p_6hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_gbl_target_residuals_nbp_vs_p_5hit = TH2F("h_gbl_target_residuals_nbp_vs_p_5hit", 
	"Non-Bend Plane Residuals at Target vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_gbl_target_residuals_nbp_vs_p_6hit = TH2F("h_gbl_target_residuals_nbp_vs_p_6hit", 
	"Non-Bend Plane Residuals at Target vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 160, -20, 20)

h_target_residuals_bp_vs_p_5hit = TH2F("h_target_residuals_bp_vs_p_5hit", 
		"Bend Plane Residuals at the Target - 5 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_target_residuals_bp_vs_p_6hit = TH2F("h_target_residuals_bp_vs_p_6hit", 
		"Bend Plane Residuals at the Target - 6 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_sp_residuals_bp_vs_p_5hit = TH2F("h_sp_residuals_bp_vs_p_5hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_sp_residuals_bp_vs_p_6hit = TH2F("h_sp_residuals_bp_vs_p_6hit", 
	"Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 240, -30, 30)
h_sp_residuals_nbp_vs_p_5hit = TH2F("h_sp_residuals_nbp_vs_p_5hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_sp_residuals_nbp_vs_p_6hit = TH2F("h_sp_residuals_nbp_vs_p_6hit", 
	"Non-Bend Plane Residuals at Scoring Plane vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_target_residuals_nbp_vs_p_5hit = TH2F("h_target_residuals_nbp_vs_p_5hit", 
	"Non-Bend Plane Residuals at Target vs Momentum - 5 Hit Tracks", 10, 0, 2.5, 160, -20, 20)
h_target_residuals_nbp_vs_p_6hit = TH2F("h_target_residuals_nbp_vs_p_6hit", 
	"Non-Bend Plane Residuals at Target vs Momentum - 6 Hit Tracks", 10, 0, 2.5, 160, -20, 20)


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

		# Get an Ecal scoring plane hit from the event
		ecal_sp_hit = hps_event.getEcalScoringPlaneHit(ecal_sp_hit_n)

		# Get the associated particle
		particle = ecal_sp_hit.getParticle()

		# Get the track associated with the particle
		track = particle.getTracks().At(0)

		# Find the GBL track that corresponds to the Ecal scoring 
		# plane hit
		for gbl_track_n in xrange(0, hps_event.getNumberOfGblTracks()):

			# Get a GBL track from the event
			gbl_track = hps_event.getGblTrack(gbl_track_n)

			# If the GBL seed track and the track associated with
			# the Ecal scoring plane hit do not match, continue
			# onto the next track
			if(gbl_track.getSeedD0() != track.getD0()): continue

			h_track_d0.Fill(gbl_track.getSeedD0())
			h_gbl_track_d0.Fill(gbl_track.getD0())

			gbl_track_parameters[0] = gbl_track.getD0()	
			gbl_track_parameters[1] = gbl_track.getPhi()
			gbl_track_parameters[2] = gbl_track.getKappa()
			gbl_track_parameters[3] = gbl_track.getZ0()
			gbl_track_parameters[4] = math.tan(math.pi/2 - gbl_track.getTheta())

			track_parameters[0] = track.getD0()
			track_parameters[1] = track.getPhi()
			track_parameters[2] = track.getOmega()
			track_parameters[3] = track.getZ0()
			track_parameters[4] = track.getTanLambda()

			extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(gbl_track_parameters)
			gbl_track_position = extrapolator.extrapolate_track(ecal_sp_hit.getPosition()[2])

			gbl_delta_x = ecal_sp_hit.getPosition()[0] - gbl_track_position[0]
			gbl_delta_y = ecal_sp_hit.getPosition()[1] - gbl_track_position[1]
			
			extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(track_parameters)
			track_position = extrapolator.extrapolate_track(ecal_sp_hit.getPosition()[2])
			
			delta_x = ecal_sp_hit.getPosition()[0] - track_position[0]
			delta_y = ecal_sp_hit.getPosition()[1] - track_position[1]
				
			momentum = particle.getMomentum()[0]*particle.getMomentum()[0]
			momentum += particle.getMomentum()[1]*particle.getMomentum()[1]
			momentum += particle.getMomentum()[2]*particle.getMomentum()[2]
			momentum = math.sqrt(momentum)

			#print "Track hits: " + str(track.getNumberOfHits())
			if track.getNumberOfHits() == 5:
				h_sp_residuals_bp_vs_p_5hit.Fill(momentum, delta_x)
				h_sp_residuals_nbp_vs_p_5hit.Fill(momentum, delta_y)
				h_gbl_sp_residuals_bp_vs_p_5hit.Fill(momentum, gbl_delta_x)
				h_gbl_sp_residuals_nbp_vs_p_5hit.Fill(momentum, gbl_delta_y)
			elif track.getNumberOfHits() == 6:
				h_sp_residuals_bp_vs_p_6hit.Fill(momentum, delta_x)
				h_sp_residuals_nbp_vs_p_6hit.Fill(momentum, delta_y)
				h_gbl_sp_residuals_bp_vs_p_6hit.Fill(momentum, gbl_delta_x)
				h_gbl_sp_residuals_nbp_vs_p_6hit.Fill(momentum, gbl_delta_y)

			extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(gbl_track_parameters)
			gbl_track_position = extrapolator.extrapolate_track(.1)

			gbl_delta_x = 0 - gbl_track_position[0]
			gbl_delta_y = 0 - gbl_track_position[1]
			
			extrapolator = HpsTrackExtrapolator.HpsTrackExtrapolator(track_parameters)
			track_position = extrapolator.extrapolate_track(.1)
			
			delta_x = 0 - track_position[0]
			delta_y = 0 - track_position[1]
				
			#print "Track hits: " + str(track.getNumberOfHits())
			if track.getNumberOfHits() == 5:
				h_target_residuals_bp_vs_p_5hit.Fill(momentum, delta_x)
				h_target_residuals_nbp_vs_p_5hit.Fill(momentum, delta_y)
				h_gbl_target_residuals_bp_vs_p_5hit.Fill(momentum, gbl_delta_x)
				h_gbl_target_residuals_nbp_vs_p_5hit.Fill(momentum, gbl_delta_y)
			elif track.getNumberOfHits() == 6:
				h_target_residuals_bp_vs_p_6hit.Fill(momentum, delta_x)
				h_target_residuals_nbp_vs_p_6hit.Fill(momentum, delta_y)
				h_gbl_target_residuals_bp_vs_p_6hit.Fill(momentum, gbl_delta_x)
				h_gbl_target_residuals_nbp_vs_p_6hit.Fill(momentum, gbl_delta_y)


canvas.Print("track_parameters.pdf[")

#--- Fit the D0 plots ---#
#------------------------#

track_d0 = RooRealVar("track_d0", "D0 (mm)", -1.5, 1.5)
track_d0_mean = RooRealVar("track_d0_mean", "Track D0 Mean", 0, -1.5, 1.5)
track_d0_sigma = RooRealVar("track_d0_sigma", "Track D0 Sigma", .5, 0, 2)
track_d0_gaussian = RooGaussModel("track_d0_gaussian", "Track D0 Gaussian", track_d0, track_d0_mean, track_d0_sigma)

track_d0_data = RooDataHist("track_d0_data", "Track D0 Data", RooArgList(track_d0), h_track_d0)
track_d0_plot = track_d0.frame()
track_d0_plot.SetTitle("")
track_d0_data.plotOn(track_d0_plot, RooFit.MarkerColor(kAzure+2), RooFit.LineColor(kAzure+2))

track_d0_gaussian.fitTo(track_d0_data)
track_d0_gaussian.plotOn(track_d0_plot, RooFit.LineColor(kAzure+2))

gbl_track_d0 = RooRealVar("gbl_track_d0", "D0 (mm)", -1.5, 1.5)
gbl_track_d0_mean = RooRealVar("gbl_track_d0_mean", "Track D0 Mean", 0, -1.5, 1.5)
gbl_track_d0_sigma = RooRealVar("gbl_track_d0_sigma", "Track D0 Sigma", .5, 0, 2)
gbl_track_d0_gaussian = RooGaussModel("gbl_track_d0_gaussian", "Track D0 Gaussian", gbl_track_d0, gbl_track_d0_mean, gbl_track_d0_sigma)

gbl_track_d0_data = RooDataHist("track_d0_data", "Track D0 Data", RooArgList(gbl_track_d0), h_gbl_track_d0)
gbl_track_d0_plot = gbl_track_d0.frame()
gbl_track_d0_plot.SetTitle("")
gbl_track_d0_data.plotOn(gbl_track_d0_plot, RooFit.MarkerColor(kRed+1), RooFit.LineColor(kRed+1))

gbl_track_d0_gaussian.fitTo(gbl_track_d0_data)
gbl_track_d0_gaussian.plotOn(gbl_track_d0_plot, RooFit.LineColor(kRed+1))

track_d0_plot.Draw()
gbl_track_d0_plot.Draw("same")
legend.Draw()
canvas.Print("track_parameters.pdf(")

canvas.Print("track_parameters.pdf]")

#--- Resolutions ---#
#-------------------#

sp_bp_residuals = RooRealVar("sp_residuals_bp", "Bend Plane Residuals (mm)", -30, 30)

bp_core_mean = RooRealVar("bp_core_mean", "bp_core_mean", 0, -2, 2)
bp_core_sigma = RooRealVar("bp_core_sigma", "bp_core_sigma", 1, 0, 5)
bp_core_gaussian = RooGaussModel("bp_core_gaussian", "bp_core_gaussian", sp_bp_residuals, bp_core_mean, bp_core_sigma)

bp_tail_mean = RooRealVar("bp_tail_mean", "bp_tail_mean", 0, -10, 10)
bp_tail_sigma = RooRealVar("bp_tail_sigma", "bp_tail_sigma", 9, 0, 15)
bp_tail_gaussian = RooGaussModel("bp_tail_gaussian", "bp_tail_gaussian", sp_bp_residuals, bp_tail_mean, bp_tail_sigma)

bp_core_fraction = RooRealVar("bp_core_fraction", "Fraction of Core", .8, 0, 1)
resolution_model = RooAddModel("resolution_model", "resolution_model", RooArgList(bp_core_gaussian, bp_tail_gaussian),
		RooArgList(bp_core_fraction))

canvas.Print("sp_residuals.pdf[")

h_sp_residuals_bp_vs_p_5hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

g_bp_pos_resolution_vs_p_5hit = TGraphErrors()

for bin in xrange(1, h_sp_residuals_bp_vs_p_5hit.GetXaxis().GetNbins()):

	projection = h_sp_residuals_bp_vs_p_5hit.ProjectionY("", bin, bin)
	
	if projection.GetEntries() < 100: continue

	sp_bp_residuals_data = RooDataHist("sp_bp_residuals_data", "Bend Plane Residuals Data", RooArgList(sp_bp_residuals), projection)
	sp_bp_residuals_plot = sp_bp_residuals.frame()
	sp_bp_residuals_data.plotOn(sp_bp_residuals_plot)
	
	resolution_model.fitTo(sp_bp_residuals_data)
	resolution_model.plotOn(sp_bp_residuals_plot)

	g_bp_pos_resolution_vs_p_5hit.SetPoint(bin-1, bin*.25, bp_core_sigma.getValV())
	g_bp_pos_resolution_vs_p_5hit.SetPointError(bin-1, 0, bp_core_sigma.getError())

	sp_bp_residuals_plot.Draw()
	canvas.Print("sp_residuals.pdf(")
	
	sp_bp_residuals_plot.remove()
	sp_bp_residuals_plot.remove()


g_bp_pos_resolution_vs_p_5hit.Draw("*Ae")
canvas.Print("sp_residuals.pdf(")

sp_nbp_residuals = RooRealVar("sp_residuals_nbp", "Non-Bend Plane Residuals (mm)", -20, 20)

nbp_core_mean = RooRealVar("nbp_core_mean", "nbp_core_mean", 0, -2, 2)
nbp_core_sigma = RooRealVar("nbp_core_sigma", "nbp_core_sigma", 1, 0, 5)
nbp_core_gaussian = RooGaussModel("nbp_core_gaussian", "nbp_core_gaussian", sp_nbp_residuals, nbp_core_mean, nbp_core_sigma)

nbp_tail_mean = RooRealVar("nbp_tail_mean", "nbp_tail_mean", 0, -10, 10)
nbp_tail_sigma = RooRealVar("nbp_tail_sigma", "nbp_tail_sigma", 9, 0, 15)
nbp_tail_gaussian = RooGaussModel("nbp_tail_gaussian", "nbp_tail_gaussian", sp_nbp_residuals, nbp_tail_mean, nbp_tail_sigma)

nbp_core_fraction = RooRealVar("nbp_core_fraction", "Fraction of Core", .8, 0, 1)
resolution_model = RooAddModel("resolution_model", "resolution_model", RooArgList(nbp_core_gaussian, nbp_tail_gaussian), RooArgList(nbp_core_fraction))

h_gbl_sp_residuals_bp_vs_p_5hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

#--- 5 hit plots ---#
#-------------------#

h_sp_residuals_nbp_vs_p_5hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

g_nbp_pos_resolution_vs_p_5hit = TGraphErrors()

for bin in xrange(1, h_sp_residuals_nbp_vs_p_5hit.GetXaxis().GetNbins()):

	projection = h_sp_residuals_nbp_vs_p_5hit.ProjectionY("", bin, bin)
	
	if projection.GetEntries() < 100: continue

	sp_nbp_residuals_data = RooDataHist("sp_nbp_residuals_data", "Non-Bend Plane Residuals Data", RooArgList(sp_nbp_residuals), projection)
	sp_nbp_residuals_plot = sp_nbp_residuals.frame()
	sp_nbp_residuals_data.plotOn(sp_nbp_residuals_plot)
	
	resolution_model.fitTo(sp_nbp_residuals_data)
	resolution_model.plotOn(sp_nbp_residuals_plot)

	g_nbp_pos_resolution_vs_p_5hit.SetPoint(bin-1, bin*.25, nbp_core_sigma.getValV())
	g_nbp_pos_resolution_vs_p_5hit.SetPointError(bin-1, 0, nbp_core_sigma.getError())

	sp_nbp_residuals_plot.Draw()
	canvas.Print("sp_residuals.pdf(")
	
	sp_nbp_residuals_plot.remove()
	sp_nbp_residuals_plot.remove()

g_nbp_pos_resolution_vs_p_5hit.Draw("A*")
canvas.Print("sp_residuals.pdf(")

h_gbl_sp_residuals_nbp_vs_p_5hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

h_sp_residuals_bp_vs_p_6hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

h_gbl_sp_residuals_bp_vs_p_6hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

h_gbl_sp_residuals_nbp_vs_p_6hit.Draw("colz")
canvas.Print("sp_residuals.pdf(")

canvas.Print("sp_residuals.pdf]")

canvas.Print("target_residuals.pdf[")

h_target_residuals_bp_vs_p_5hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_target_residuals_nbp_vs_p_5hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_gbl_target_residuals_bp_vs_p_5hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_gbl_target_residuals_nbp_vs_p_5hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_target_residuals_bp_vs_p_6hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_target_residuals_nbp_vs_p_6hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_gbl_target_residuals_bp_vs_p_6hit.Draw("colz")
canvas.Print("target_residuals.pdf(")
h_gbl_target_residuals_nbp_vs_p_6hit.Draw("colz")

canvas.Print("target_residuals.pdf]")

