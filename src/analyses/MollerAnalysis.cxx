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
      class_name("MollerAnalysis") {

}

MollerAnalysis::~MollerAnalysis() {
    delete plotter;
}

void MollerAnalysis::initialize() { 
    this->bookHistograms(); 
}

void MollerAnalysis::processEvent(HpsEvent* event) { 

    // Only look at pair1 triggers
    if (!event->isPair1Trigger()) return;

    // Only look at events that have two clusters
    if (event->getNumberOfEcalClusters() != 2) return;

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

    if (abs(delta_cluster_time) > 3) return;

    plotter->get2DHistogram("first v second cluster energy - N clusters == 2, time cut")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("first v second cluster time - N clusters == 2, time cut")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());


    // Only look at events that have two tracks 
    if (event->getNumberOfTracks() >= 2) return;

    plotter->get1DHistogram("cluster energy - N tracks >= 2")->Fill(first_cluster->getEnergy());
    plotter->get1DHistogram("cluster energy - N tracks >= 2")->Fill(second_cluster->getEnergy());

    SvtTrack* first_track = NULL;
    SvtTrack* second_track = NULL;
    std::cout << "Matching tracks" << std::endl; 
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
        if (this->isMatch(first_cluster, event->getTrack(track_n))) { 
            first_track = event->getTrack(track_n);
            std::cout << "First track match found" << std::endl;
        } else if (this->isMatch(second_cluster, event->getTrack(track_n))) { 
            second_track = event->getTrack(track_n);
            std::cout << "Second track match found" << std::endl;
        }
    }
    
    if (first_track == NULL || second_track == NULL) return;
    
    plotter->get2DHistogram("first v second cluster energy - matched")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    plotter->get2DHistogram("first v second cluster time - matched")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());
    
    /*
    matcher->findAllMatches(event); 
    SvtTrack* first_track = matcher->getMatchingTrack(first_cluster);
    SvtTrack* second_track = matcher->getMatchingTrack(second_cluster);



    double p0 = AnalysisUtils::getMagnitude(first_track->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(second_track->getMomentum()); 

    plotter->get2DHistogram("p[e-] v p[e-] - matched")->Fill(p0, p1);

    // Require that both tracks are negatively charged
    if ((first_track->getCharge() + second_track->getCharge()) != -2) return;
    */

    /*
    SvtTrack* first_track = event->getTrack(0);
    SvtTrack* second_track = event->getTrack(1); 



    if (first_track->isTopTrack()) { 
        plotter->get1DHistogram("p top")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom")->Fill(p0);
    }

    if (second_track->isTopTrack()) { 
        plotter->get1DHistogram("p top")->Fill(p0);
    } else { 
        plotter->get1DHistogram("p bottom")->Fill(p0);
    }


    plotter->get1DHistogram("delta Track Time - all")->Fill(first_track->getTrackTime() - second_track->getTrackTime());

    plotter->get1DHistogram("delta cos(theta)")->Fill(TrackExtrapolator::getCosTheta(first_track) - 
            TrackExtrapolator::getCosTheta(second_track));

    // Require that the tracks are in opposite volumes
    if ((first_track->isTopTrack() && second_track->isTopTrack()) 
            || (first_track->isBottomTrack() && second_track->isBottomTrack())) return;

    // Require that both tracks are matched to a cluster
    if (matcher->getMatchingCluster(first_track) == NULL || 
            matcher->getMatchingCluster(second_track) == NULL) return;
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->Fill(p0, p1);
    plotter->get1DHistogram("delta Track Time - matched")->Fill(first_track->getTrackTime() - second_track->getTrackTime());

    first_cluster = matcher->getMatchingCluster(first_track);
    second_cluster = matcher->getMatchingCluster(second_track);

    plotter->get2DHistogram("first v second cluster time")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    plotter->get2DHistogram("first v second cluster energy - matched")->Fill(
            first_cluster->getEnergy(), 
            second_cluster->getEnergy());

    delta_cluster_time = first_cluster->getClusterTime() 
        - second_cluster->getClusterTime();

    plotter->get1DHistogram("delta clusters - matched")->Fill(delta_cluster_time);

    if (abs(delta_cluster_time) > 10) return;

    plotter->get2DHistogram("first v second cluster time - pass time cut")->Fill(
            first_cluster->getClusterTime(), 
            second_cluster->getClusterTime());

    // Require that the sum of the momentum be above and below 
    // some reasonable value
    //if ((p0+p1) < .80*1.056) return;
    //if ((p0+p1) > 1.2*1.056) return;

    //if (abs(first_track->getTrackTime() - second_track->getTrackTime()) > 4) return;

    //if (TrackExtrapolator::getCosTheta(first_track) - 
    //        TrackExtrapolator::getCosTheta(second_track) < 0.04) return;

    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->Fill(p0, p1);
    plotter->get1DHistogram("delta Track Time - pass cuts")->Fill(first_track->getTrackTime() - second_track->getTrackTime());
    plotter->get1DHistogram("delta clusters - pass cuts")->Fill(delta_cluster_time); */

}

void MollerAnalysis::finalize() { 
    plotter->saveToPdf("moller_analysis.pdf");
}


void MollerAnalysis::bookHistograms() {

    plotter->build1DHistogram("p top", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");
    plotter->build1DHistogram("p bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("p (GeV)");

    plotter->build1DHistogram("cluster energy - N clusters == 2", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster Energy (GeV)");
    plotter->build1DHistogram("cluster energy - N tracks >= 2", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster Energy (GeV)");

    plotter->build1DHistogram("delta Track Time - all", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");
    plotter->build1DHistogram("delta Track Time - matched", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");
    plotter->build1DHistogram("delta Track Time - pass cuts", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Track time [ns]");

    plotter->build1DHistogram("delta clusters - matched", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Cluster time [ns]");

    plotter->build1DHistogram("delta clusters - pass cuts", 100, -10, 10)->GetXaxis()->SetTitle("#Delta Cluster time [ns]");

    plotter->build1DHistogram("delta cos(theta)", 40, -0.1, 0.1); 

    plotter->build2DHistogram("p[e-] v p[e-] - matched", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - matched")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

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

    plotter->build2DHistogram("p[e-] v p[e-] - Pass Cuts", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetXaxis()->SetTitle("p[e-] [GeV]"); 
    plotter->get2DHistogram("p[e-] v p[e-] - Pass Cuts")->GetYaxis()->SetTitle("p[e-] [GeV]"); 

    plotter->build2DHistogram("cluster x v extrapolated track x - top", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y - top", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - top", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - top", 50, -25, 25);
    plotter->build2DHistogram("cluster x v extrapolated track x - bottom", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y - bottom", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x - bottom", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y - bottom", 50, -25, 25);
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
    } else {
        plotter->get2DHistogram("cluster x v extrapolated track x - bottom")->Fill(cluster_pos[0], 
                track_pos_at_cluster_shower_max[0]);
        plotter->get2DHistogram("cluster y v extrapolated track y - bottom")->Fill(cluster_pos[1], 
                track_pos_at_cluster_shower_max[1]);

        plotter->get1DHistogram("cluster x - extrapolated track x - bottom")->Fill(cluster_pos[0] 
                - track_pos_at_cluster_shower_max[0]);
        plotter->get1DHistogram("cluster y - extrapolated track y - bottom")->Fill(cluster_pos[1] 
                - track_pos_at_cluster_shower_max[1]);
    }
    
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (abs(cluster_pos[1] - track_pos_at_cluster_shower_max[1]) > 30) return false;
    
    /*if (cluster_pos[0] - track_pos_at_cluster_shower_max[0] > 30 ||
            cluster_pos[0] - track_pos_at_cluster_shower_max[0] < -30) return false;
    
    if (cluster_pos[1] - track_pos_at_cluster_shower_max[1] > 30
            || cluster_pos[1] - track_pos_at_cluster_shower_max[1] < -30) return false;*/

    return true;
}
