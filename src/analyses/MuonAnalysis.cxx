/**
 * @file MuonAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date July 10, 2015
 */

#include <MuonAnalysis.h>

MuonAnalysis::MuonAnalysis()
    : plotter(new Plotter()),
      class_name("MuonAnalysis") {

}

MuonAnalysis::~MuonAnalysis() {
    delete plotter;
}

void MuonAnalysis::initialize() { 
    this->bookHistograms(); 
}

void MuonAnalysis::processEvent(HpsEvent* event) { 


    // Only look at pair1 triggers
    //if (!event->isPair1Trigger()) return;
    // Only look at single1 triggers
    if (!event->isSingle1Trigger()) return;

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = AnalysisUtils::getClusterPair(event);
    if (pair.size() != 2) return;
    
    plotter->get2DHistogram("cluster pair energy")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);

    plotter->get2DHistogram("cluster y v cluster x")->Fill(pair[1]->getPosition()[0], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster energy sum")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

    // Require that both clusters are small
    if (pair[0]->getEcalHits()->GetEntriesFast() >= 3 
            || pair[1]->getEcalHits()->GetEntriesFast() >= 3) return;

    plotter->get2DHistogram("cluster pair energy - cuts: size")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt - cuts: size")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time - cuts: size")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position - cuts: size")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: size")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);
    plotter->get2DHistogram("cluster y v cluster x - cuts: size")->Fill(pair[0]->getPosition()[0], 
            pair[0]->getPosition()[1]);
    plotter->get2DHistogram("cluster y v cluster x - cuts: size")->Fill(pair[1]->getPosition()[0], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster energy sum - cuts: size")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

    // Require that the two clusters are in opposite volumes.  This cut should
    // eventually become part of the standard pair requirement.
    if (pair[0]->getPosition()[1]*pair[1]->getPosition()[1] > 0) return;

    // Require that one cluster is beams left and the other beams right
    if (pair[0]->getPosition()[0]*pair[1]->getPosition()[0] > 0) return;

    plotter->get2DHistogram("cluster pair energy - cuts: size, fiducial")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt - cuts: size, fiducial")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time - cuts: size, fiducial")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position - cuts: size, fiducial")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: size, fiducial")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);
    plotter->get2DHistogram("cluster y v cluster x - cuts: size, fiducial")->Fill(pair[1]->getPosition()[0], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster energy sum - cuts: size, fiducial")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

    SvtTrack* first_track = NULL;
    SvtTrack* second_track = NULL; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        if (AnalysisUtils::isMatch(pair[0], event->getTrack(track_n))) { 
            first_track = event->getTrack(track_n);
        } else if (AnalysisUtils::isMatch(pair[1], event->getTrack(track_n))) { 
            second_track = event->getTrack(track_n);
        }
    }

    if (first_track == NULL && second_track == NULL) {
        return;
    } else if (first_track == NULL || second_track == NULL) { 
        return;
    }
    
    plotter->get2DHistogram("cluster pair energy - matched")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt - matched")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time - matched")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position - matched")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - matched")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);
    plotter->get2DHistogram("cluster y v cluster x - matched")->Fill(pair[1]->getPosition()[0], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster energy sum - matched")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());

}

void MuonAnalysis::finalize() {

    plotter->saveToPdf("muon_analysis.pdf");
    plotter->saveToRootFile("muon_analysis.root");
}


void MuonAnalysis::bookHistograms() {

    //---------------------//
    //-- Cluster energy ---//
    //---------------------//

    plotter->build2DHistogram("cluster pair energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: size", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: size")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: size")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: size, fiducial", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: size, fiducial")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build1DHistogram("cluster energy sum", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy sum (GeV)");
    plotter->build1DHistogram("cluster energy sum - cuts: size", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy sum (GeV)");
    plotter->build1DHistogram("cluster energy sum - cuts: size, fiducial", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy sum (GeV)");
    plotter->build1DHistogram("cluster energy sum - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy sum (GeV)");

    // Cluster time //

    plotter->build1DHistogram("cluster pair time dt", 40, -10, 10);

    plotter->build2DHistogram("cluster pair time", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster pair time dt - cuts: size", 40, -10, 10);
    
    plotter->build2DHistogram("cluster pair time - cuts: size", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - cuts: size")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: size")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster pair time dt - cuts: size, fiducial", 40, -10, 10);

    plotter->build2DHistogram("cluster pair time - cuts: size, fiducial", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - cuts: size, fiducial")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster pair time dt - matched", 40, -10, 10);

    plotter->build2DHistogram("cluster pair time - matched", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - matched")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - matched")->GetYaxis()->SetTitle("Cluster time (ns)");
    
    //-------------------------//
    //--- Cluster positions ---//
    //-------------------------//

    // All clusters
    plotter->build2DHistogram("cluster x position", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x position")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y position", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y position")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y position")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster y v cluster x", 200, -200, 200, 100, -100, 100);

    // Cluster passing size cuts
    plotter->build2DHistogram("cluster x position - cuts: size", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - cuts: size")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x position - cuts: size")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y position - cuts: size", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y position - cuts: size")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y position - cuts: size")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster y v cluster x - cuts: size", 200, -200, 200, 100, -100, 100);

    // Cluster passing fiducial cuts
    plotter->build2DHistogram("cluster x position - cuts: size, fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x position - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y position - cuts: size, fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y position - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y position - cuts: size, fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster y v cluster x - cuts: size, fiducial", 200, -200, 200, 100, -100, 100);

    plotter->build2DHistogram("cluster x position - matched", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - matched")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x position - matched")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y position - matched", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y position - matched")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y position - matched")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster y v cluster x - matched", 200, -200, 200, 100, -100, 100);


} 

std::string MuonAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
