/**
 * @file MollerAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date July 10, 2015
 */

#include <MollerAnalysis.h>

MollerAnalysis::MollerAnalysis()
    : plotter(new Plotter()),
      matcher(new TrackClusterMatcher()),
      class_name("MollerAnalysis"), 
      total_events(0),
      total_pair_trigger_events(0),  
      total_pair_events(0), 
      total_two_cluster_events(0) {

}

MollerAnalysis::~MollerAnalysis() {
    delete plotter;
}

void MollerAnalysis::initialize() { 
    this->bookHistograms(); 
}

void MollerAnalysis::processEvent(HpsEvent* event) { 

    total_events++;

    // Only look at pair1 triggers
    if (!event->isPair1Trigger()) return;
    total_pair_trigger_events++;

    this->getClusterPair(event);

    // Only look at events that have two clusters
    if (event->getNumberOfEcalClusters() != 2) return;
    total_two_cluster_events++; 

    EcalCluster* first_cluster = event->getEcalCluster(0);
    EcalCluster* second_cluster = event->getEcalCluster(1);

    plotter->get2DHistogram("first v second cluster energy - N clusters == 2")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());


    plotter->get2DHistogram("first v second cluster time - N clusters == 2")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    plotter->get1DHistogram("cluster energy - N clusters == 2")->Fill(first_cluster->getEnergy());
    plotter->get1DHistogram("cluster energy - N clusters == 2")->Fill(second_cluster->getEnergy());

    double delta_cluster_time = first_cluster->getClusterTime() - second_cluster->getClusterTime();

    if (std::abs(delta_cluster_time) > 2.5) return;

    plotter->get2DHistogram("first v second cluster energy - N clusters == 2, time cut")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("first v second cluster time - N clusters == 2, time cut")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    // Only look at events that have two tracks 
    if (event->getNumberOfTracks() < 2) return;

    plotter->get1DHistogram("cluster energy - N tracks >= 2")->Fill(first_cluster->getEnergy());
    plotter->get1DHistogram("cluster energy - N tracks >= 2")->Fill(second_cluster->getEnergy());

    //std::cout << "Matching " << event->getNumberOfTracks() << " tracks" << std::endl;
    SvtTrack* first_track = NULL;
    SvtTrack* second_track = NULL; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        //std::cout << "Trying to find match for track " << track_n << std::endl;
        if (this->isMatch(first_cluster, event->getTrack(track_n))) { 
            first_track = event->getTrack(track_n);
            //std::cout << "First track match found" << std::endl;
        } else if (this->isMatch(second_cluster, event->getTrack(track_n))) { 
            second_track = event->getTrack(track_n);
            //std::cout << "Second track match found" << std::endl;
        }
    }
    
    if (first_track == NULL && second_track == NULL) {

        for (int first_track_n = 0; first_track_n < event->getNumberOfTracks(); ++first_track_n) { 
            first_track = event->getTrack(first_track_n);
            for (int second_track_n = 0; second_track_n < event->getNumberOfTracks(); ++second_track_n) { 
                if (first_track_n == second_track_n) continue;
                second_track = event->getTrack(second_track_n);

                plotter->get2DHistogram("first v second track time - no match")->Fill(first_track->getTrackTime(), second_track->getTrackTime()); 
                plotter->get2DHistogram("first v second theta - no match")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
                        std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

                double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
                double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

                plotter->get2DHistogram("p v theta - no match")->Fill(p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
                plotter->get2DHistogram("p v theta - no match")->Fill(p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));
                plotter->get2DHistogram("p[e] v p[e] - no match")->Fill(p0, p1);
            } 
        }
        //std::cout << "One of the two tracks is NULL" << std::endl;
        return;
    }

    if (first_track == NULL || second_track == NULL) { 
        for (int first_track_n = 0; first_track_n < event->getNumberOfTracks(); ++first_track_n) { 
            first_track = event->getTrack(first_track_n);
            for (int second_track_n = 0; second_track_n < event->getNumberOfTracks(); ++second_track_n) { 
                if (first_track_n == second_track_n) continue;
                second_track = event->getTrack(second_track_n);

                plotter->get2DHistogram("first v second track time - single match")->Fill(first_track->getTrackTime(), second_track->getTrackTime()); 
                plotter->get2DHistogram("first v second theta - single match")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
                        std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

                double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
                double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

                plotter->get2DHistogram("p v theta - single match")->Fill(p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
                plotter->get2DHistogram("p v theta - single match")->Fill(p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));
                plotter->get2DHistogram("p[e] v p[e] - single match")->Fill(p0, p1);
            } 
        }
        return;        
    }

    plotter->get2DHistogram("first v second cluster energy - matched")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("first v second cluster time - matched")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

    plotter->get1DHistogram("delta Track Time - matched")->Fill(first_track->getTrackTime() - second_track->getTrackTime());
    plotter->get2DHistogram("first v second track time - matched")->Fill(first_track->getTrackTime(), second_track->getTrackTime());

    /*std::cout << "cos(theta): " << TrackExtrapolator::getCosTheta(first_track) 
              << " theta: " << 3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track)) 
              << " pi/2 - theta " << 3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track)) << std::endl;*/

    plotter->get2DHistogram("first v second theta - matched")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));


    plotter->get2DHistogram("p v theta - matched")->Fill(p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched")->Fill(p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    if (std::abs(first_track->getTrackTime() - second_track->getTrackTime()) > 5) return;

    plotter->get2DHistogram("p[e] v p[e] - matched")->Fill(p0, p1);

    if (first_track->isTopTrack()) { 
        plotter->get1DHistogram("p top - matched")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom - matched")->Fill(p0);
    }

    if (second_track->isTopTrack()) { 
        plotter->get1DHistogram("p top - matched")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom - matched")->Fill(p0);
    }
   
    plotter->get2DHistogram("first v second cluster x position - matched")->Fill(first_cluster->getPosition()[0], 
            second_cluster->getPosition()[0]);
    plotter->get2DHistogram("first v second cluster y position - matched")->Fill(first_cluster->getPosition()[1], 
            second_cluster->getPosition()[1]);

    // Require that both tracks are negatively charged
    if ((first_track->getCharge() + second_track->getCharge()) != -2) return;

    if (p0 > .85 || p1 > .85) return;

    plotter->get2DHistogram("p[e-] v p[e-] - matched")->Fill(p0, p1);
    plotter->get2DHistogram("first v second track time - matched, e-e-")->Fill(first_track->getTrackTime(), second_track->getTrackTime());

    plotter->get2DHistogram("first v second theta - matched, e-e-")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("p v theta - matched, e-e-")->Fill(p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched, e-e-")->Fill(p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("first v second cluster energy - matched, e-e-")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("first v second cluster x position - matched, e-e-")->Fill(first_cluster->getPosition()[0], 
            second_cluster->getPosition()[0]);
    plotter->get2DHistogram("first v second cluster y position - matched, e-e-")->Fill(first_cluster->getPosition()[1], 
            second_cluster->getPosition()[1]);

}

void MollerAnalysis::finalize() {

    //TF1* p_v_theta = new TF1("p_v_theta", "sqrt(2*(.000510/1.056)*(1.056/x -1 + .000510/x))", 0, 1.2);

    std::cout << "Total number of events: " << total_events << std::endl;
    std::cout << "Pair trigger events: " << total_pair_trigger_events << " / " << total_events 
              << " = " << (total_pair_trigger_events/total_events)*100 << " % " << std::endl;  
    std::cout << "Pair events: " << total_pair_events << " / " << total_events 
              << " = " << (total_pair_events/total_events)*100 << " % " << std::endl;  
    std::cout << "Two cluster events: " << total_two_cluster_events << " / " << total_events 
              << " = " << (total_two_cluster_events/total_events)*100 << " % " << std::endl;  


    plotter->saveToPdf("moller_analysis.pdf");
    plotter->saveToRootFile("moller_analysis.root");
}


void MollerAnalysis::bookHistograms() {

    plotter->build1DHistogram("p top - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build1DHistogram("cluster energy - N clusters == 2", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->build1DHistogram("cluster energy - N tracks >= 2", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->build1DHistogram("cluster energy - N tracks == 1", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build1DHistogram("delta Track Time - matched", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");

    plotter->build1DHistogram("delta theta - matched", 40, -0.1, 0.1); 
    plotter->build1DHistogram("delta theta - matched, e+e-", 40, -0.1, 0.1); 

    plotter->build2DHistogram("first v second track time - no match", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("first v second track time - no match")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("first v second track time - no match")->GetYaxis()->SetTitle("Track time [ns]");

    plotter->build2DHistogram("first v second track time - single match", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("first v second track time - single match")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("first v second track time - single match")->GetYaxis()->SetTitle("Track time [ns]");

    plotter->build2DHistogram("first v second track time - matched", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("first v second track time - matched")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("first v second track time - matched")->GetYaxis()->SetTitle("Track time [ns]");
    
    plotter->build2DHistogram("first v second theta - matched", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("first v second theta - matched")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("first v second theta - matched")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("first v second theta - no match", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("first v second theta - no match")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("first v second theta - no match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("first v second theta - single match", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("first v second theta - single match")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("first v second theta - single match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - no match", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - no match")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - no match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - single match", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - single match")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - single match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p[e] v p[e] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - matched")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - matched")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e] v p[e] - no match", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - no match")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - no match")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e] v p[e] - single match", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - single match")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - single match")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("first v second track time - matched, e-e-", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("first v second track time - matched, e-e-")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("first v second track time - matched, e-e-")->GetYaxis()->SetTitle("Track time [ns]");

    plotter->build2DHistogram("p v theta - matched, e-e-", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched, e-e-")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched, e-e-")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("first v second theta - matched, e-e-", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("first v second theta - matched, e-e-")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("first v second theta - matched, e-e-")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("first v second cluster time - N clusters == 2", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("first v second cluster time - N clusters == 2")->GetXaxis()->SetTitle("First cluster time (ns)");
    plotter->get2DHistogram("first v second cluster time - N clusters == 2")->GetYaxis()->SetTitle("Second cluster time (ns)");

    plotter->build2DHistogram("first v second cluster time - N clusters == 2, time cut", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("first v second cluster time - N clusters == 2, time cut")->GetXaxis()->SetTitle("First cluster time (ns)");
    plotter->get2DHistogram("first v second cluster time - N clusters == 2, time cut")->GetYaxis()->SetTitle("Second cluster time (ns)");
    
    plotter->build2DHistogram("first v second cluster time - matched", 160, 20, 100, 160, 20, 100);
    plotter->get2DHistogram("first v second cluster time - matched")->GetXaxis()->SetTitle("First cluster time (ns)");
    plotter->get2DHistogram("first v second cluster time - matched")->GetYaxis()->SetTitle("Second cluster time (ns)");

    plotter->build2DHistogram("first v second cluster energy - N clusters == 2", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("first v second cluster energy - N clusters == 2")->GetYaxis()->SetTitle("First cluster energy (GeV)");
    plotter->get2DHistogram("first v second cluster energy - N clusters == 2")->GetYaxis()->SetTitle("Second cluster energy (GeV)");

    plotter->build2DHistogram("first v second cluster energy - N clusters == 2, time cut", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("first v second cluster energy - N clusters == 2, time cut")->GetYaxis()->SetTitle("First cluster energy (GeV)");
    plotter->get2DHistogram("first v second cluster energy - N clusters == 2, time cut")->GetYaxis()->SetTitle("Second cluster energy (GeV)");

    plotter->build2DHistogram("first v second cluster energy - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("first v second cluster energy - matched")->GetYaxis()->SetTitle("First cluster energy (GeV)");
    plotter->get2DHistogram("first v second cluster energy - matched")->GetYaxis()->SetTitle("Second cluster energy (GeV)");

    plotter->build2DHistogram("first v second cluster x position - matched", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("first v second cluster x position - matched")->GetYaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("first v second cluster x position - matched")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("first v second cluster y position - matched", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("first v second cluster y position - matched")->GetYaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("first v second cluster y position - matched")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("first v second cluster energy - matched, e-e-", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("first v second cluster energy - matched, e-e-")->GetYaxis()->SetTitle("First cluster energy (GeV)");
    plotter->get2DHistogram("first v second cluster energy - matched, e-e-")->GetYaxis()->SetTitle("Second cluster energy (GeV)");

    plotter->build2DHistogram("first v second cluster x position - matched, e-e-", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("first v second cluster x position - matched, e-e-")->GetYaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("first v second cluster x position - matched, e-e-")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("first v second cluster y position - matched, e-e-", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("first v second cluster y position - matched, e-e-")->GetYaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("first v second cluster y position - matched, e-e-")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster x v extrapolated track x - top", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster x v extrapolated track x - top - electrons", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster x v extrapolated track x - top - positrons", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("p v extrapolated track x - top", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("p v cluster x - top", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("cluster energy v cluster x - top", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("cluster x - track x v e/p - top", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster x v e/p - top", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v e/p - top", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v extrapolated track y - top", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - top", 100, -100, 100);
    plotter->build1DHistogram("cluster y - extrapolated track y - top", 100, -100, 100);
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom - electrons", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom - positrons", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("p v extrapolated track x - bottom", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("p v cluster x - bottom", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("cluster energy v cluster x - bottom", 50, 0, 1.5, 100, -100, 100);
    plotter->build2DHistogram("cluster x - track x v e/p - bottom", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster x v e/p - bottom", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v e/p - bottom", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - bottom", 100, -100, 100);
    plotter->build1DHistogram("cluster y - extrapolated track y - bottom", 100, -100, 100);

    plotter->build2DHistogram("cluster x v extrapolated track x - top - matched", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y - top - matched", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - top - matched", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - top - matched", 50, -25, 25);
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom - matched", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom - matched", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - bottom - matched", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - bottom - matched", 50, -25, 25);
} 

std::string MollerAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}


bool MollerAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();
    /*std::cout << "[ TagProbeAnalysis ]: ECal cluster position: " 
        << " x: " << cluster_pos[0] 
        << " y: " << cluster_pos[1] 
        << " z: " << cluster_pos[2]
        << std::endl;*/

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_cluster_shower_max 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);
    /*std::cout << "[ TagProbeAnalysis ]: Track position at shower max: " 
         << " x: " << track_pos_at_cluster_shower_max[0] 
         << " y: " << track_pos_at_cluster_shower_max[1] 
         << " z: " << track_pos_at_cluster_shower_max[2]
         << std::endl;*/ 


    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    // If the track and cluster are in opposite volumes, then they can't 
    // be a match
    if (cluster_pos[1]*track_pos_at_cluster_shower_max[1] < 0) return false;

    if (track->isTopTrack()) { 
        plotter->get2DHistogram("cluster x v extrapolated track x - top")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - top")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - top")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);

        plotter->get2DHistogram("p v extrapolated track x - top")->Fill(p, track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("p v cluster x - top")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster energy v cluster x - top")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - top")->Fill(cluster_pos[0] - track_pos_at_cluster_shower_max[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - top")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - top")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 

        if (track->getCharge() < 0) { 
            plotter->get2DHistogram("cluster x v extrapolated track x - top - electrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        } else {
            plotter->get2DHistogram("cluster x v extrapolated track x - top - positrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        }
    
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
        
        plotter->get2DHistogram("p v extrapolated track x - bottom")->Fill(p, track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("p v cluster x - bottom")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster energy v cluster x - bottom")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - bottom")->Fill(cluster_pos[0] - track_pos_at_cluster_shower_max[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - bottom")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - bottom")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 
        
        if (track->getCharge() < 0) { 
            plotter->get2DHistogram("cluster x v extrapolated track x - bottom - electrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        } else {
            plotter->get2DHistogram("cluster x v extrapolated track x - bottom - positrons")->Fill(cluster_pos[0], 
                    track_pos_at_cluster_shower_max[0]);
        }
    }
    
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (std::abs(cluster_pos[0] - track_pos_at_cluster_shower_max[0]) > 20) return false;

    if (cluster->getEnergy()/p < .5) return false;

    if (track->isTopTrack()) { 
        plotter->get2DHistogram("cluster x v extrapolated track x - top - matched")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top - matched")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - top - matched")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - top - matched")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom - matched")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom - matched")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom - matched")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom - matched")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
    }

    /*if (cluster_pos[0] - track_pos_at_cluster_shower_max[0] > 30 ||
            cluster_pos[0] - track_pos_at_cluster_shower_max[0] < -30) return false;
    
    if (cluster_pos[1] - track_pos_at_cluster_shower_max[1] > 30
            || cluster_pos[1] - track_pos_at_cluster_shower_max[1] < -30) return false;*/

    //std::cout << "Track and cluster are a match" << std::endl;
    return true;
}

std::vector<EcalCluster*> MollerAnalysis::getClusterPair(HpsEvent* event) { 

    EcalCluster* first_cluster;
    EcalCluster* second_cluster;

    // Loop through all clusters in an event
    for (int first_cluster_n = 0; first_cluster_n < event->getNumberOfEcalClusters(); ++first_cluster_n) { 
    
        // Get an Ecal cluster from the event
        first_cluster = event->getEcalCluster(first_cluster_n);

        // Loop through the rest of the clusters and make pairs
        for (int second_cluster_n = (first_cluster_n + 1); second_cluster_n < event->getNumberOfEcalClusters();
                ++second_cluster_n) { 
            
            // Get another Ecal cluster from the event
            second_cluster = event->getEcalCluster(second_cluster_n);

            // Check if the two clusters can be considered a 'good pair'. This is done by requiring 
            // the clusters to satisfy a series of cuts
             
            double delta_cluster_time = first_cluster->getClusterTime() - second_cluster->getClusterTime();
            if (std::abs(delta_cluster_time) > 2.5) continue;

            total_pair_events++;
            break;
        }
    }

    std::vector<EcalCluster*> pair;
    pair.push_back(first_cluster);
    pair.push_back(second_cluster);

    return pair; 
}
