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
    //if (!event->isPair1Trigger()) return;
    // Only look at single1 triggers
    if (!event->isSingle1Trigger()) return;
    total_pair_trigger_events++;

    // Get a "good" pair from the event.  If a good pair isn't found, skip
    // the event.
    std::vector<EcalCluster*> pair = AnalysisUtils::getClusterPair(event);
    if (pair.size() != 2) return;
    total_pair_events++;
    
    plotter->get2DHistogram("cluster pair energy")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);

    // Require that the two clusters are in opposite volumes.  This cut should
    // eventually become part of the standard pair requirement.
    if (pair[0]->getPosition()[1]*pair[1]->getPosition()[1] > 0) return;

    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->Fill(pair[0]->getEnergy(), pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair time dt - cuts: fiducial")->Fill(pair[0]->getClusterTime() - pair[1]->getClusterTime());
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->Fill(pair[0]->getClusterTime(), pair[1]->getClusterTime()); 
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);


    SvtTrack* first_track = NULL;
    SvtTrack* second_track = NULL; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        //std::cout << "Trying to find match for track " << track_n << std::endl;
        if (this->isMatch(pair[0], event->getTrack(track_n))) { 
            first_track = event->getTrack(track_n);
            //std::cout << "First track match found" << std::endl;
        } else if (this->isMatch(pair[1], event->getTrack(track_n))) { 
            second_track = event->getTrack(track_n);
            //std::cout << "Second track match found" << std::endl;
        }
    }

    if (first_track == NULL && second_track == NULL) {
        return;
    } else if (first_track == NULL || second_track == NULL) { 
        return;
    }

    plotter->get2DHistogram("cluster pair energy - matched")->Fill(
            pair[0]->getEnergy(), 
            pair[1]->getEnergy());

    plotter->get2DHistogram("cluster pair time - matched")->Fill(
            pair[0]->getClusterTime(), 
            pair[1]->getClusterTime());

    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

    plotter->get1DHistogram("delta Track Time - matched")->Fill(
            first_track->getTrackTime() - second_track->getTrackTime());
    plotter->get2DHistogram("track time - matched")->Fill(first_track->getTrackTime(), second_track->getTrackTime());

    plotter->get2DHistogram("track theta - matched")->Fill(
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("p v theta - matched")->Fill(
            p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched")->Fill(
            p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));


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
   
    plotter->get2DHistogram("cluster x position - matched")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - matched")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster pair dx - matched")->Fill(pair[0]->getPosition()[0] - pair[1]->getPosition()[0]);

    // Require that both tracks are negatively charged
    if ((first_track->getCharge() + second_track->getCharge()) != -2) return;

    plotter->get2DHistogram("p[e-] v p[e-] - matched")->Fill(p0, p1);
    plotter->get2DHistogram("track time - matched, e-e-")->Fill(first_track->getTrackTime(),
            second_track->getTrackTime());

    plotter->get2DHistogram("track theta - matched, e-e-")->Fill(
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("p v theta - matched, e-e-")->Fill(p0, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched, e-e-")->Fill(p1, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("cluster pair energy - matched, e-e-")->Fill(
            pair[0]->getEnergy(), 
            pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - matched, e-e-")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - matched, e-e-")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster pair dx - matched, e-e-")->Fill(pair[0]->getPosition()[0] - pair[1]->getPosition()[0]);

    if (std::abs(first_track->getTrackTime() - second_track->getTrackTime()) > 4) return;

    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->Fill(p0, p1);
    
    plotter->get2DHistogram("track theta - matched, e-e-, time")->Fill(
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("p v theta - matched, e-e-, time")->Fill(p0, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched, e-e-, time")->Fill(p1, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));
    
    plotter->get2DHistogram("track time - matched, e-e-, time")->Fill(first_track->getTrackTime(),
            second_track->getTrackTime());

   
    //if (p0 > .85 || p1 > .85) return;
    if (p0+p1 > 1.1*1.056) return;
    
    if (p0+p1 < .85*1.056) return;


    plotter->get2DHistogram("p[e-] v p[e-] - matched, time, p")->Fill(p0, p1);
    
    plotter->get2DHistogram("track theta - matched, e-e-, time, p")->Fill(
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("p v theta - matched, e-e-, time, p")->Fill(p0, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched, e-e-, time, p")->Fill(p1, 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

    plotter->get2DHistogram("cluster pair energy - moller")->Fill(
            pair[0]->getEnergy(), 
            pair[1]->getEnergy());

    plotter->get1DHistogram("cluster pair energy sum - moller")->Fill(pair[0]->getEnergy() + pair[1]->getEnergy());
    plotter->get1DHistogram("cluster pair energy diff - moller")->Fill(pair[0]->getEnergy() - pair[1]->getEnergy());

    plotter->get2DHistogram("cluster x position - moller")->Fill(pair[0]->getPosition()[0], 
            pair[1]->getPosition()[0]);
    plotter->get2DHistogram("cluster y position - moller")->Fill(pair[0]->getPosition()[1], 
            pair[1]->getPosition()[1]);

    plotter->get2DHistogram("cluster position - moller")->Fill(pair[0]->getPosition()[0], 
            pair[0]->getPosition()[1]);
    plotter->get2DHistogram("cluster position - moller")->Fill(pair[1]->getPosition()[0], 
            pair[1]->getPosition()[1]);

    plotter->get1DHistogram("cluster pair dx - moller")->Fill(pair[0]->getPosition()[0] - pair[1]->getPosition()[0]);

    //plotter->get1DHistogram("angle between clusters - moller")->Fill(atan(

    //double p_expected = sqrt(2*(.000510/1.056)*(1.056/x -1 + .000510/x))
    //double p_residual = 

    //plotter->get1DHistogram("p residual")->Fill();

    //--- Two cluster events ---//
    //--------------------------//

    // Only look at events that have two clusters
    /*if (event->getNumberOfEcalClusters() != 2) return;
    total_two_cluster_events++; 

    EcalCluster* first_cluster = event->getEcalCluster(0);
    EcalCluster* second_cluster = event->getEcalCluster(1);

    plotter->get2DHistogram("cluster pair energy - N clusters == 2")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());


    plotter->get2DHistogram("cluster pair time - N clusters == 2")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    plotter->get1DHistogram("cluster pair energy - N clusters == 2")->Fill(first_cluster->getEnergy());
    plotter->get1DHistogram("cluster pair energy - N clusters == 2")->Fill(second_cluster->getEnergy());

    double delta_cluster_time = first_cluster->getClusterTime() - second_cluster->getClusterTime();

    if (std::abs(delta_cluster_time) > 2.5) return;

    plotter->get2DHistogram("cluster pair energy - N clusters == 2, time cut")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("cluster pair time - N clusters == 2, time cut")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    // Only look at events that have two tracks 
    if (event->getNumberOfTracks() < 2) return;

    plotter->get1DHistogram("cluster pair energy - N tracks >= 2")->Fill(first_cluster->getEnergy());
    plotter->get1DHistogram("cluster pair energy - N tracks >= 2")->Fill(second_cluster->getEnergy());

    //std::cout << "Matching " << event->getNumberOfTracks() << " tracks" << std::endl;
    first_track = NULL;
    second_track = NULL; 
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

                plotter->get2DHistogram("track time - no match")->Fill(first_track->getTrackTime(), second_track->getTrackTime()); 
                plotter->get2DHistogram("track theta - no match")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
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

                plotter->get2DHistogram("track time - single match")->Fill(first_track->getTrackTime(), second_track->getTrackTime()); 
                plotter->get2DHistogram("track theta - single match")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
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
*/
    /*
    plotter->get2DHistogram("cluster pair energy - matched")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("cluster pair time - matched")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

    plotter->get1DHistogram("delta Track Time - matched")->Fill(first_track->getTrackTime() - second_track->getTrackTime());
    plotter->get2DHistogram("track time - matched")->Fill(first_track->getTrackTime(), second_track->getTrackTime());

    plotter->get2DHistogram("track theta - matched")->Fill(std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))), 
            std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));


    plotter->get2DHistogram("p v theta - matched")->Fill(p0, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(first_track))));
    plotter->get2DHistogram("p v theta - matched")->Fill(p1, std::abs(3.14159/2 - acos(TrackExtrapolator::getCosTheta(second_track))));

*/
    
    



}

void MollerAnalysis::finalize() {

    //TF1* p_v_theta = new TF1("p_v_theta", "sqrt(2*(.000510/1.056)*(1.056/x -1 + .000510/x))", 0, 1.2);

    std::cout << "Total number of events: " << total_events << std::endl;
    std::cout << "Total trigger events: " << total_pair_trigger_events << " / " << total_events 
              << " = " << (total_pair_trigger_events/total_events)*100 << " % " << std::endl;  
    std::cout << "Total events with a pair: " << total_pair_events << " / " << total_events 
              << " = " << (total_pair_events/total_events)*100 << " % " << std::endl;  
    std::cout << "Two cluster events: " << total_two_cluster_events << " / " << total_events 
              << " = " << (total_two_cluster_events/total_events)*100 << " % " << std::endl;  


    plotter->saveToPdf("moller_analysis.pdf");
    plotter->saveToRootFile("moller_analysis.root");
}


void MollerAnalysis::bookHistograms() {

    //---------------------//
    //-- Cluster energy ---//
    //---------------------//

    plotter->build2DHistogram("cluster pair energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetYaxis()->SetTitle("Cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("First cluster pair energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched")->GetYaxis()->SetTitle("Second cluster pair energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - matched, e-e-", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - matched, e-e-")->GetYaxis()->SetTitle("First cluster pair energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - matched, e-e-")->GetYaxis()->SetTitle("Second cluster pair energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - moller", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - moller")->GetYaxis()->SetTitle("First cluster pair energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - moller")->GetYaxis()->SetTitle("Second cluster pair energy (GeV)");

    plotter->build1DHistogram("cluster pair energy sum - moller", 50, 0, 1.5);
    plotter->build1DHistogram("cluster pair energy diff - moller", 50, -1, 1);

    plotter->build2DHistogram("cluster pair energy - N clusters == 2", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - N clusters == 2")->GetYaxis()->SetTitle(
            "First cluster pair energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - N clusters == 2")->GetYaxis()->SetTitle(
            "Second cluster pair energy (GeV)");

    // Cluster time //

    plotter->build1DHistogram("cluster pair time dt", 40, -10, 10);

    plotter->build2DHistogram("cluster pair time", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster pair time dt - cuts: fiducial", 40, -10, 10);

    plotter->build2DHistogram("cluster pair time - cuts: fiducial", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - n clusters == 2", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - n clusters == 2")->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - n clusters == 2")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - matched", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - matched")->GetXaxis()->SetTitle("First cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - matched")->GetYaxis()->SetTitle("Second cluster time (ns)");

    // Track momentum

    plotter->build1DHistogram("p top - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build2DHistogram("p[e] v p[e] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - matched")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - matched")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e] v p[e] - single match", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - single match")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - single match")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e] v p[e] - no match", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e] v p[e] - no match")->GetXaxis()->SetTitle("p[e] [GeV]"); 
    plotter->get2DHistogram("p[e] v p[e] - no match")->GetYaxis()->SetTitle("p[e] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched, time", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched, time, p", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time, p")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched, time, p")->GetYaxis()->SetTitle("p[e-] [GeV]");

    plotter->build1DHistogram("p residual", 50, -40, 40); 

    // Track time

    plotter->build1DHistogram("delta Track Time - matched", 100, -10, 10)->GetXaxis()
        ->SetTitle("#Delta Track time [ns]");

    plotter->build2DHistogram("track time - matched", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("track time - matched")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("track time - matched")->GetYaxis()->SetTitle("Track time [ns]");

    plotter->build2DHistogram("track time - matched, e-e-", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("track time - matched, e-e-")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("track time - matched, e-e-")->GetYaxis()->SetTitle("Track time [ns]");

    plotter->build2DHistogram("track time - matched, e-e-, time", 100, -10, 10, 100, -10, 10);
    plotter->get2DHistogram("track time - matched, e-e-, time")->GetXaxis()->SetTitle("Track time [ns]");
    plotter->get2DHistogram("track time - matched, e-e-, time")->GetYaxis()->SetTitle("Track time [ns]");

    // Track theta

    plotter->build1DHistogram("delta theta - matched", 40, -0.1, 0.1); 
    
    plotter->build1DHistogram("delta theta - matched, e+e-", 40, -0.1, 0.1); 

    plotter->build2DHistogram("track theta - matched", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - matched, e-e-", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched, e-e-")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched, e-e-")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - matched, e-e-, time", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched, e-e-, time")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched, e-e-, time")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - matched, e-e-, time, p", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - matched, e-e-, time, p")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - matched, e-e-, time, p")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - single match", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - single match")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - single match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("track theta - no match", 40, 0, 0.1, 40, 0, 0.1);
    plotter->get2DHistogram("track theta - no match")->GetXaxis()->SetTitle("#theta");
    plotter->get2DHistogram("track theta - no match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched, e-e-", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched, e-e-")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched, e-e-")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched, e-e-, time", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched, e-e-, time")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched, e-e-, time")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - matched, e-e-, time, p", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - matched, e-e-, time, p")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - matched, e-e-, time, p")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - single match", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - single match")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - single match")->GetYaxis()->SetTitle("#theta");

    plotter->build2DHistogram("p v theta - no match", 50, 0, 1.5, 40, 0, 0.1);
    plotter->get2DHistogram("p v theta - no match")->GetXaxis()->SetTitle("p (GeV)");
    plotter->get2DHistogram("p v theta - no match")->GetYaxis()->SetTitle("#theta");

    
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

    // Cluster passing fiducial cuts
    plotter->build2DHistogram("cluster x position - cuts: fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x position - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y position - cuts: fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y position - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    // Clusters matched to tracks
    plotter->build2DHistogram("cluster x position - matched", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - matched")->GetYaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - matched")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("cluster y position - matched", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - matched")->GetYaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - matched")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build1DHistogram("cluster pair dx - matched", 50, -150, 150);

    plotter->build2DHistogram("cluster x position - matched, e-e-", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - matched, e-e-")->GetYaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - matched, e-e-")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("cluster y position - matched, e-e-", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - matched, e-e-")->GetYaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - matched, e-e-")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build1DHistogram("cluster pair dx - matched, e-e-", 50, -150, 150);

    plotter->build2DHistogram("cluster x position - moller", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x position - moller")->GetYaxis()->SetTitle("First cluster x position (mm)");
    plotter->get2DHistogram("cluster x position - moller")->GetYaxis()->SetTitle("Second cluster x position (mm)");

    plotter->build2DHistogram("cluster y position - moller", 100, -100, 100, 100, -100, 100);
    plotter->get2DHistogram("cluster y position - moller")->GetYaxis()->SetTitle("First cluster y position (mm)");
    plotter->get2DHistogram("cluster y position - moller")->GetYaxis()->SetTitle("Second cluster y position (mm)");

    plotter->build2DHistogram("cluster position - moller", 200, -200, 200, 100, -100, 100);

    plotter->build1DHistogram("cluster pair dx - moller", 50, -150, 150);

    //------------------------------------//
    //--- Extrapolated track positions ---//
    //------------------------------------//

    // All tracks and clusters
    
    // Top
    plotter->build2DHistogram("cluster x v extrapolated track x - top", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("p v extrapolated track x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("p v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster pair energy v cluster x - top", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster x - track x v e/p - top", 200, -200, 200, 40, 0, 2);
    plotter->build2DHistogram("cluster x v e/p - top", 200, -200, 200, 40, 0, 2);
    plotter->build1DHistogram("cluster x - extrapolated track x - top", 200, -200, 200);
   
     plotter->build1DHistogram("cluster y - extrapolated track y - top", 100, -100, 100);
    plotter->build2DHistogram("cluster y v e/p - top", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v extrapolated track y - top", 100, -100, 100, 100, -100, 100);

    // Bottom
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("p v extrapolated track x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("p v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster pair energy v cluster x - bottom", 50, 0, 1.5, 200, -200, 200);
    plotter->build2DHistogram("cluster x - track x v e/p - bottom", 200, -200, 200, 40, 0, 2);
    plotter->build2DHistogram("cluster x v e/p - bottom", 200, -200, 200, 40, 0, 2);
    plotter->build1DHistogram("cluster x - extrapolated track x - bottom", 200, -200, 200);

    plotter->build2DHistogram("cluster y v e/p - bottom", 100, -100, 100, 40, 0, 2);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster y - extrapolated track y - bottom", 100, -100, 100);

    // Matched tracks

    // Top
    plotter->build2DHistogram("cluster x v extrapolated track x - top - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("cluster y v extrapolated track y - top - matched", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - top - matched", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - top - matched", 50, -25, 25);
    
    // Bottom
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom - matched", 200, -200, 200, 200, -200, 200);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom - matched", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - bottom - matched", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - bottom - matched", 50, -25, 25);
} 

std::string MollerAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}


bool MollerAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

    // Check that the track and cluster are in the same detector volume.
    // If not, thre is no way they can match.
    if (track->isTopTrack() && cluster->getPosition()[1] < 0
            || track->isBottomTrack() && cluster->getPosition()[1] > 0) return false;
    
    // Get the cluster position
    std::vector<double> cluster_pos = cluster->getPosition();

    // Extrapolate the track to the Ecal cluster position 
    std::vector<double> track_pos_at_ecal 
        = TrackExtrapolator::extrapolateTrack(track, cluster_pos[2]);

    double p = AnalysisUtils::getMagnitude(track->getMomentum());

    double delta_x = cluster_pos[0] - track_pos_at_ecal[0];
    double delta_y = cluster_pos[1] - track_pos_at_ecal[1];

    if (track->isTopTrack()) { 
        plotter->get2DHistogram("cluster x v extrapolated track x - top")->Fill(cluster_pos[0], 
                track_pos_at_ecal[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top")->Fill(cluster_pos[1], 
                track_pos_at_ecal[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - top")->Fill(delta_x); 
        plotter->get1DHistogram("cluster y - extrapolated track y - top")->Fill(delta_y); 

        plotter->get2DHistogram("p v extrapolated track x - top")->Fill(p, track_pos_at_ecal[0]);
        plotter->get2DHistogram("p v cluster x - top")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - top")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - top")->Fill(delta_x, cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - top")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - top")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 
    
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom")->Fill(cluster_pos[0], 
                track_pos_at_ecal[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom")->Fill(cluster_pos[1], 
                track_pos_at_ecal[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom")->Fill(delta_x);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom")->Fill(delta_y);
        
        plotter->get2DHistogram("p v extrapolated track x - bottom")->Fill(p, track_pos_at_ecal[0]);
        plotter->get2DHistogram("p v cluster x - bottom")->Fill(p, cluster_pos[0]);
        plotter->get2DHistogram("cluster pair energy v cluster x - bottom")->Fill(cluster->getEnergy(), cluster_pos[0]);
        plotter->get2DHistogram("cluster x - track x v e/p - bottom")->Fill(delta_x,
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster x v e/p - bottom")->Fill(cluster_pos[0],
               cluster->getEnergy()/p); 
        plotter->get2DHistogram("cluster y v e/p - bottom")->Fill(cluster_pos[1],
               cluster->getEnergy()/p); 
        
    }
    
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if ((track->isTopTrack() && (delta_x > 14 || delta_x < -18)) ||
        (track->isBottomTrack() && (delta_x > 9 || delta_x < -21))) return false;

    if ((track->isTopTrack() && (delta_y > 14 || delta_y < -14)) ||
        (track->isBottomTrack() && (delta_y > 14 || delta_y < -14))) return false;

    //if (cluster->getEnergy()/p < .5) return false;

    if (track->isTopTrack()) { 
        plotter->get2DHistogram("cluster x v extrapolated track x - top - matched")->Fill(cluster_pos[0], 
                track_pos_at_ecal[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - top - matched")->Fill(cluster_pos[1], 
                track_pos_at_ecal[1]);
        plotter->get1DHistogram("cluster x - extrapolated track x - top - matched")->Fill(delta_x);
        plotter->get1DHistogram("cluster y - extrapolated track y - top - matched")->Fill(delta_y);
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom - matched")->Fill(cluster_pos[0], 
                track_pos_at_ecal[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom - matched")->Fill(cluster_pos[1], 
                track_pos_at_ecal[1]);
        plotter->get1DHistogram("cluster x - extrapolated track x - bottom - matched")->Fill(delta_x);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom - matched")->Fill(delta_y);
    }

    return true;
}

