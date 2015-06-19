
#include <SimpleTrackingEfficiencyAnalysis.h>

SimpleTrackingEfficiencyAnalysis::SimpleTrackingEfficiencyAnalysis() 
    : track(NULL),
      cluster(NULL),
      plotter(new Plotter()),
      cluster_energy_low_threshold(.7 /* GeV */),
      cluster_energy_high_threshold(1.15 /* GeV */),
      cuts_enabled(true),
      class_name("SimpleTrackingEfficiencyAnalysis") {

}

SimpleTrackingEfficiencyAnalysis::~SimpleTrackingEfficiencyAnalysis() {
}

void SimpleTrackingEfficiencyAnalysis::initialize() { 
    this->bookHistograms();
}

void SimpleTrackingEfficiencyAnalysis::processEvent(HpsEvent* event) { 
  
    // Only look at single 1 triggers
    if (!event->isSingle1Trigger()) return;

    // Loop over all clusters in an event and try to match a track to them
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
  
        // Get the cluster from the event 
        cluster = event->getEcalCluster(cluster_n);
    
        // Reset the cut flags
        bool pass_time_cut = false;
        bool pass_cluster_size_cut = false;
        bool pass_energy_cut = false;

        // Fill the cluster information for all events
        double cluster_energy = cluster->getEnergy();
        double cluster_time = cluster->getClusterTime();
        plotter->get1DHistogram("cluster energy")->Fill(cluster_energy);
        plotter->get1DHistogram("cluster time")->Fill(cluster_time);
        plotter->get2DHistogram("cluster energy v cluster time - all")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size - all")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());

        // Get the seed hit of the cluster
        EcalHit* seed_hit = cluster->getSeed();

        // Make the same plots for the top and bottom Ecal volumes 
        if (seed_hit->getYCrystalIndex() > 0) { 
            plotter->get1DHistogram("cluster energy - top")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - top")->Fill(cluster_time);
        } else { 
            plotter->get1DHistogram("cluster energy - bottom")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - bottom")->Fill(cluster_time);
        }  

        plotter->get2DHistogram("cluster count - all")->Fill(
                seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);

        if (!isEdgeCrystal(seed_hit)) { 
            plotter->get1DHistogram("cluster energy - no edge")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - no edge")->Fill(cluster_time);
        }
 
        // Check that the cluster passes the time requirement
        if (passClusterTimeCut(cluster)) {
            plotter->get2DHistogram("cluster energy v cluster time - pass time cut")->Fill(cluster_energy, cluster_time);
            pass_time_cut = true;
        }

        // 
        if (passClusterSizeCut(cluster)) { 
            pass_cluster_size_cut = true; 
        } 


        // Check that the cluster passes the energy requirement
        if (passEnergyCut(cluster)) {
           pass_energy_cut = true; 
        }

        if (pass_time_cut && pass_cluster_size_cut && pass_energy_cut) { 
            plotter->get2DHistogram("cluster energy v cluster time - fee")->Fill(cluster_energy, cluster_time);
            plotter->get2DHistogram("cluster energy v cluster size - fee")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
            plotter->get2DHistogram("cluster count - fee")->Fill(
                    seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);
        }

        bool match = false;
        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
            
            // Get a track from the event
            track = event->getTrack(track_n);

            plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
            plotter->get1DHistogram("cluster time - track time")->Fill(
                    cluster->getClusterTime() - track->getTrackTime());


            // Calculate the momentum magnitude and transverse momentum
            std::vector<double> p = track->getMomentum();
            double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

            plotter->get1DHistogram("E/p")->Fill(cluster_energy/p_mag);
        
            if (!isEdgeCrystal(seed_hit)) { 
                plotter->get1DHistogram("E/p - no edge")->Fill(cluster_energy/p_mag);
            }


            if (isMatch(cluster, track)) {
                plotter->get2DHistogram("tracking efficiency - all")->Fill(seed_hit->getXCrystalIndex(), 
                    seed_hit->getYCrystalIndex(), 1);
                plotter->get1DHistogram("tracking efficiency - cluster energy")->Fill(cluster->getEnergy(), 1);
                plotter->get1DHistogram("tracking efficiency - cluster time")->Fill(cluster->getClusterTime(), 1);
                
                plotter->get1DHistogram("doca - match")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - match")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - match")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - match")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - match")->Fill(track->getTanLambda());

                if (track->getCharge() < 0) { 
                    plotter->get1DHistogram("doca - electron - match")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - electron - match")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - electron - match")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - electron - match")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - electron - match")->Fill(track->getTanLambda());
                } else { 
                    plotter->get1DHistogram("doca - positron - match")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - positron - match")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - positron - match")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - positron - match")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - positron - match")->Fill(track->getTanLambda());
                }

                if (pass_time_cut && pass_cluster_size_cut && pass_energy_cut) { 
                    plotter->get2DHistogram("tracking efficiency - fee")->Fill(seed_hit->getXCrystalIndex(), 
                        seed_hit->getYCrystalIndex(), 1);
                    plotter->get1DHistogram("E/p - fee")->Fill(cluster_energy/p_mag);
                }

                if (!isEdgeCrystal(seed_hit)) { 
                    plotter->get1DHistogram("tracking efficiency - cluster energy - no edge")->Fill(cluster->getEnergy(), 1);
                    plotter->get1DHistogram("tracking efficiency - cluster time - no edge")->Fill(cluster->getClusterTime(), 1);
                }
                match = true;
                break; 
            } else { 
                plotter->get1DHistogram("doca - no match")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - no match")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - no match")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - no match")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - no match")->Fill(track->getTanLambda());

                if (track->getCharge() < 0) { 
                    plotter->get1DHistogram("doca - electron - no match")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - electron - no match")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - electron - no match")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - electron - no match")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - electron - no match")->Fill(track->getTanLambda());
                } else { 
                    plotter->get1DHistogram("doca - positron - no match")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - positron - no match")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - positron - no match")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - positron - no match")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - positron - no match")->Fill(track->getTanLambda());
                }
            } 
        }
        if (pass_time_cut && pass_cluster_size_cut && pass_energy_cut && !match) { 
            plotter->get2DHistogram("cluster count - no match")->Fill(
                    seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);
        }
    }
}

void SimpleTrackingEfficiencyAnalysis::finalize() { 
    
    plotter->get2DHistogram("tracking efficiency - fee")->Divide(
            plotter->get2DHistogram("cluster count - fee"));
    plotter->get2DHistogram("tracking efficiency - all")->Divide(
            plotter->get2DHistogram("cluster count - all"));
    plotter->get1DHistogram("tracking efficiency - cluster energy")->Divide(
            plotter->get1DHistogram("cluster energy"));
    plotter->get1DHistogram("tracking efficiency - cluster time")->Divide(
            plotter->get1DHistogram("cluster time"));
    plotter->get1DHistogram("tracking efficiency - cluster energy - no edge")->Divide(
            plotter->get1DHistogram("cluster energy - no edge"));
    plotter->get1DHistogram("tracking efficiency - cluster time - no edge")->Divide(
            plotter->get1DHistogram("cluster time - no edge"));


    plotter->saveToPdf("simple_tracking_efficiency.pdf");
    plotter->saveToRootFile("simple_tracking_efficiency.root");
}

void SimpleTrackingEfficiencyAnalysis::bookHistograms() { 


    plotter->build1DHistogram("cluster energy", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time", 160, 0, 80);
    plotter->build1DHistogram("track time", 40, -20, 20);
    plotter->build1DHistogram("cluster time - track time", 100, -20, 80);
    plotter->build1DHistogram("E/p", 40, 0, 1.5);
    plotter->build2DHistogram("cluster x v extrapolated track x", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y", 50, -25, 25);
    plotter->build2DHistogram("cluster count - all", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("cluster energy v cluster size - all", 50, 0, 1.5, 10, 0, 10);
    plotter->build2DHistogram("cluster energy v cluster time - all", 50, 0, 1.5, 160, 0, 80);
    plotter->build2DHistogram("tracking efficiency - all", 47, -23, 24, 12, -6, 6);

    plotter->build1DHistogram("tracking efficiency - cluster energy", 50, 0, 1.5);
    plotter->build1DHistogram("tracking efficiency - cluster time", 160, 0, 80);

    plotter->build1DHistogram("cluster energy - top", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - top", 160, 0, 80);

    plotter->build1DHistogram("cluster energy - bottom", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - bottom", 160, 0, 80);

    plotter->build1DHistogram("cluster energy - no edge", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - no edge", 160, 0, 80);
    plotter->build1DHistogram("E/p - no edge", 40, 0, 1.5);
    plotter->build1DHistogram("tracking efficiency - cluster energy - no edge", 50, 0, 1.5);
    plotter->build1DHistogram("tracking efficiency - cluster time - no edge", 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster time - pass time cut", 50, 0, 1.5, 160, 0, 80);
    
    plotter->build2DHistogram("cluster energy v cluster time - fee", 50, 0, 1.5, 160, 0, 80);
    plotter->build2DHistogram("cluster energy v cluster size - fee", 50, 0, 1.5, 10, 0, 10);
    plotter->build2DHistogram("cluster count - fee", 47, -23, 24, 12, -6, 6);
    plotter->build2DHistogram("tracking efficiency - fee", 47, -23, 24, 12, -6, 6);
    plotter->build1DHistogram("E/p - fee", 40, 0, 1.5);

    plotter->build2DHistogram("cluster count - no match", 47, -23, 24, 12, -6, 6);
    plotter->build1DHistogram("doca - no match", 80, -10, 10);
    plotter->build1DHistogram("z0 - no match", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - no match", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - no match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - no match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - no match", 40, -0.1, 0.1);
    plotter->build1DHistogram("doca - electron - no match", 80, -10, 10);
    plotter->build1DHistogram("z0 - electron - no match", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - electron - no match", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - electron - no match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - electron - no match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - electron - no match", 40, -0.1, 0.1);
    plotter->build1DHistogram("doca - positron - no match", 80, -10, 10);
    plotter->build1DHistogram("z0 - positron - no match", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - positron - no match", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - positron - no match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - positron - no match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - positron - no match", 40, -0.1, 0.1);

    plotter->build1DHistogram("doca - match", 80, -10, 10);
    plotter->build1DHistogram("z0 - match", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - match", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - match", 40, -0.1, 0.1);
    plotter->build1DHistogram("doca - electron - match", 80, -10, 10);
    plotter->build1DHistogram("sin(phi0) - electron - match", 40, -0.2, 0.2);
    plotter->build1DHistogram("z0 - electron - match", 80, -2, 2);
    plotter->build1DHistogram("curvature - electron - match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - electron - match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - electron - match", 40, -0.1, 0.1);
    plotter->build1DHistogram("doca - positron - match", 80, -10, 10);
    plotter->build1DHistogram("z0 - positron - match", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - positron - match", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - positron - match", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - positron - match", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - positron - match", 40, -0.1, 0.1);

}

std::string SimpleTrackingEfficiencyAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}


bool SimpleTrackingEfficiencyAnalysis::passEnergyCut(EcalCluster* cluster) { 
    if (cluster->getEnergy() < cluster_energy_low_threshold
            || cluster->getEnergy() > cluster_energy_high_threshold) return false;

    return true;
}

bool SimpleTrackingEfficiencyAnalysis::passClusterTimeCut(EcalCluster* cluster) {   
    if (cluster->getClusterTime() < 43 || cluster->getClusterTime() > 49) return false;

    return true;   
}

bool SimpleTrackingEfficiencyAnalysis::passClusterSizeCut(EcalCluster* cluster) { 
    if (cluster->getEcalHits()->GetEntriesFast() < 3) return false;
    return true;
}

bool SimpleTrackingEfficiencyAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

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

    plotter->get2DHistogram("cluster x v extrapolated track x")->Fill(cluster_pos[0], 
            track_pos_at_cluster_shower_max[0]);
    plotter->get2DHistogram("cluster y v extrapolated track y")->Fill(cluster_pos[1], 
            track_pos_at_cluster_shower_max[1]);

    plotter->get1DHistogram("cluster x - extrapolated track x")->Fill(cluster_pos[0] 
            - track_pos_at_cluster_shower_max[0]);
    plotter->get1DHistogram("cluster y - extrapolated track y")->Fill(cluster_pos[0] 
            - track_pos_at_cluster_shower_max[0]);
  
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if (cluster_pos[0] - track_pos_at_cluster_shower_max[0] > 12 ||
            cluster_pos[0] - track_pos_at_cluster_shower_max[0] < -18) return false;
    
    if (cluster_pos[1] - track_pos_at_cluster_shower_max[1] > 12
            || cluster_pos[1] - track_pos_at_cluster_shower_max[1] < -18) return false;

    return true;
}

bool SimpleTrackingEfficiencyAnalysis::isEdgeCrystal(EcalHit* hit) { 
    
    int x_crystal_index = hit->getXCrystalIndex();

    if (abs(x_crystal_index) == 1) return true;

    return false;

}
