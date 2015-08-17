
#include <TrackClusterMatchingEfficiencyAnalysis.h>

TrackClusterMatchingEfficiencyAnalysis::TrackClusterMatchingEfficiencyAnalysis() 
    : track(NULL),
      cluster(NULL),
      plotter(new Plotter()),
      cluster_energy_low_threshold(.7 /* GeV */),
      cluster_energy_high_threshold(1.15 /* GeV */),
      cuts_enabled(true),
      class_name("TrackClusterMatchingEfficiencyAnalysis"),
      total_events(0),
      total_single1_triggers(0),
      total_events_with_tracks(0) {

}

TrackClusterMatchingEfficiencyAnalysis::~TrackClusterMatchingEfficiencyAnalysis() {
}

void TrackClusterMatchingEfficiencyAnalysis::initialize() { 
    this->bookHistograms();
}

void TrackClusterMatchingEfficiencyAnalysis::processEvent(HpsEvent* event) { 

    // Increment the total events counter
    total_events++; 

    // Only look at single 1 triggers
    if (!event->isSingle1Trigger()) return;
    
    // Increment the singles1 trigger counter
    total_single1_triggers++;

    if (event->getNumberOfTracks() != 0) total_events_with_tracks++;

    // Loop over all clusters in an event and try to match a track to them
    matched_tracks.clear();
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
  
        // Get the cluster from the event 
        cluster = event->getEcalCluster(cluster_n);
    
        // Reset the cut flags
        bool pass_time_cut = false;
        bool pass_cluster_size_cut = false;
        bool edge_crystal_cut = false;
        bool pass_energy_cut = false;

        // Fill the cluster information for all events
        double cluster_energy = cluster->getEnergy();
        double cluster_time = cluster->getClusterTime();
        plotter->get1DHistogram("cluster energy")->Fill(cluster_energy);
        plotter->get1DHistogram("cluster time")->Fill(cluster_time);
        plotter->get2DHistogram("cluster energy v cluster time")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size")->Fill(cluster_energy, 
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

        plotter->get2DHistogram("cluster count")->Fill(
                seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);

        plotter->get2DHistogram("cluster position")->Fill(
                cluster->getPosition()[0], cluster->getPosition()[1]);

        plotter->get2DHistogram("cluster energy v cluster y")->Fill(cluster_energy, cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster energy v crystal index - y")->Fill(cluster_energy, seed_hit->getYCrystalIndex());


        if (!isEdgeCrystal(seed_hit)) { 
            plotter->get1DHistogram("cluster energy - no edge")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - no edge")->Fill(cluster_time);
        }
 
        // Check that the cluster passes the time requirement
        if (passClusterTimeCut(cluster)) {
            pass_time_cut = true;
        }

        // Check that the cluster passes the size requirement
        if (passClusterSizeCut(cluster)) {
            pass_cluster_size_cut = true; 
        }

        // Check that the cluster seed is not at the edge of the Ecal
        if (!isEdgeCrystal(seed_hit)) { 
            edge_crystal_cut = true; 
        } 

        // Check that the cluster passes the energy requirement
        if (passEnergyCut(cluster)) {
           pass_energy_cut = true; 
        }

        if (pass_time_cut) { 
            plotter->get2DHistogram("cluster energy v cluster time - time")->Fill(cluster_energy, cluster_time);
            plotter->get2DHistogram("cluster energy v cluster y - time")->Fill(cluster_energy, cluster->getPosition()[1]);
            plotter->get2DHistogram("cluster energy v crystal index - y - time")->Fill(cluster_energy, seed_hit->getYCrystalIndex());
            plotter->get2DHistogram("cluster energy v cluster size - time")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
        } 
        
        if (pass_time_cut && pass_cluster_size_cut) { 
            plotter->get2DHistogram("cluster energy v cluster time - time, size")->Fill(cluster_energy, cluster_time);
            plotter->get2DHistogram("cluster energy v cluster y - time, size")->Fill(cluster_energy, cluster->getPosition()[1]);
            plotter->get2DHistogram("cluster energy v crystal index - y - time, size")->Fill(cluster_energy, seed_hit->getYCrystalIndex());
            plotter->get2DHistogram("cluster energy v cluster size - time, size")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
        }

        if (pass_time_cut && pass_cluster_size_cut && edge_crystal_cut) { 
            plotter->get2DHistogram("cluster energy v cluster time - time, size, no edge")->Fill(cluster_energy, cluster_time);
            plotter->get2DHistogram("cluster energy v cluster y - time, size, no edge")->Fill(cluster_energy, cluster->getPosition()[1]);
            plotter->get2DHistogram("cluster energy v crystal index - y - time, size, no edge")->Fill(
                    cluster_energy, seed_hit->getYCrystalIndex());
            plotter->get2DHistogram("cluster energy v cluster size - time, size, no edge")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
        
        } 
        
        if (pass_time_cut && pass_cluster_size_cut && pass_energy_cut) { 
            plotter->get2DHistogram("cluster energy v cluster time - fee")->Fill(cluster_energy, cluster_time);
            plotter->get2DHistogram("cluster energy v cluster size - fee")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
            plotter->get2DHistogram("cluster count - fee")->Fill(
                    seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);

            if (seed_hit->getYCrystalIndex() > 0) { 
                plotter->get1DHistogram("cluster energy - top - fee")->Fill(cluster_energy);

            } else { 
                plotter->get1DHistogram("cluster energy - bottom - fee")->Fill(cluster_energy);
            } 

            if (edge_crystal_cut) { 
                int crystal_index = (int) abs(seed_hit->getYCrystalIndex());
                if (seed_hit->getYCrystalIndex() > 0) { 
                    plotter->get1DHistogram("cluster energy - top - no edge - fee")->Fill(cluster_energy);
                    plotter->get1DHistogram("cluster count - top - fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                } else { 
                    plotter->get1DHistogram("cluster energy - bottom - no edge - fee")->Fill(cluster_energy);
                    plotter->get1DHistogram("cluster count - bottom - fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                }
            } 
        }

        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

            // Get a track from the event
            track = event->getTrack(track_n);

            if (track->isTopTrack()) { 
                plotter->get1DHistogram("doca - top")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - top")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - top")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - top")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - top")->Fill(track->getTanLambda());
            } else { 
                plotter->get1DHistogram("doca - bottom")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - bottom")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - bottom")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - bottom")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - bottom")->Fill(track->getTanLambda());
            }
        }

        bool match = false;
        for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 
            
            // Get a track from the event
            track = event->getTrack(track_n);

            plotter->get1DHistogram("track time")->Fill(track->getTrackTime());
            plotter->get1DHistogram("cluster time - track time")->Fill(
                    cluster->getClusterTime() - track->getTrackTime());

            // Calculate the momentum magnitude 
            std::vector<double> p = track->getMomentum();
            double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

            plotter->get1DHistogram("E/p")->Fill(cluster_energy/p_mag);
        
            if (!isEdgeCrystal(seed_hit)) { 
                plotter->get1DHistogram("E/p - no edge")->Fill(cluster_energy/p_mag);
            }

            // Check if the track has already been matched to a cluster.  If it has, 
            // skip the track.  Otherwise, try and match the track to the cluster.
            if (std::find(matched_tracks.begin(), matched_tracks.end(), track) == matched_tracks.end() 
                    && isMatch(cluster, track)) {

                matched_tracks.push_back(track);

                plotter->get2DHistogram("cluster count - matched")->Fill(
                        seed_hit->getXCrystalIndex(), 
                        seed_hit->getYCrystalIndex(), 1);
                
                plotter->get1DHistogram("cluster energy - matched")->Fill(cluster->getEnergy(), 1);
                plotter->get1DHistogram("cluster time - matched")->Fill(cluster->getClusterTime(), 1);

                plotter->get1DHistogram("E/p - matched - fee")->Fill(cluster_energy/p_mag);
               
                if (seed_hit->getYCrystalIndex() > 0) { 
                    plotter->get1DHistogram("cluster energy - matched - top")->Fill(cluster->getEnergy(), 1);
                    plotter->get1DHistogram("cluster time - matched - top")->Fill(cluster->getClusterTime(), 1);
                } else { 
                    plotter->get1DHistogram("cluster energy - matched - bottom")->Fill(cluster->getEnergy(), 1);
                    plotter->get1DHistogram("cluster time - matched - bottom")->Fill(cluster->getClusterTime(), 1);
                }  

                if (track->isTopTrack()) { 
                    plotter->get1DHistogram("doca - top - pass cuts")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - top - pass cuts")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - top - pass cuts")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - top - pass cuts")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - top - pass cuts")->Fill(track->getTanLambda());
                } else { 
                    plotter->get1DHistogram("doca - bottom - pass cuts")->Fill(track->getD0());
                    plotter->get1DHistogram("z0 - bottom - pass cuts")->Fill(track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - bottom - pass cuts")->Fill(sin(track->getPhi0()));
                    plotter->get1DHistogram("curvature - bottom - pass cuts")->Fill(track->getOmega());
                    plotter->get1DHistogram("tan_lambda - bottom - pass cuts")->Fill(track->getTanLambda());
                }

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
                    plotter->get2DHistogram("cluster count - matched - fee")->Fill(
                            seed_hit->getXCrystalIndex(), 
                            seed_hit->getYCrystalIndex(), 1);
                    plotter->get1DHistogram("E/p - matched - fee")->Fill(cluster_energy/p_mag);
               
                    if (seed_hit->getYCrystalIndex() > 0) { 
                        plotter->get1DHistogram("cluster energy - matched - top - fee")->Fill(cluster_energy);
                    } else { 
                        plotter->get1DHistogram("cluster energy - matched - bottom - fee")->Fill(cluster_energy);
                    } 

                    if (edge_crystal_cut) { 
                        int crystal_index = (int) abs(seed_hit->getYCrystalIndex());
                        if (seed_hit->getYCrystalIndex() > 0) { 
                            plotter->get1DHistogram("cluster energy - matched - top - no edge - fee")->Fill(cluster_energy);
                            plotter->get1DHistogram("cluster count - matched - top - fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                        } else { 
                            plotter->get1DHistogram("cluster energy - matched - bottom - no edge - fee")->Fill(cluster_energy);
                            plotter->get1DHistogram("cluster count - matched - bottom - fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                        } 
                    }
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

void TrackClusterMatchingEfficiencyAnalysis::finalize() { 

    std::cout << "//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//" << std::endl;  
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of events: " << total_events << std::endl;
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of singles1 triggers: " 
              << total_single1_triggers << std::endl;
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of events with tracks: "
              << total_events_with_tracks << std::endl;
    std::cout << "//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//" << std::endl;  

    plotter->get2DHistogram("cluster count - matched")->Divide(
            plotter->get2DHistogram("cluster count"));
    plotter->get2DHistogram("cluster count - matched - fee")->Divide(
            plotter->get2DHistogram("cluster count - fee"));

    // Cluster energy 
    plotter->setGraphType("asymm");
    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency"))->Divide(
        plotter->get1DHistogram("cluster energy - matched"),
        plotter->get1DHistogram("cluster energy"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - top"),
        plotter->get1DHistogram("cluster energy - top"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - bottom"),
        plotter->get1DHistogram("cluster energy - bottom"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - top - fee"),
        plotter->get1DHistogram("cluster energy - top - fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - bottom - fee"),
        plotter->get1DHistogram("cluster energy - bottom - fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - no edge - fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - top - no edge - fee"),
        plotter->get1DHistogram("cluster energy - top - no edge - fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - no edge - fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - bottom - no edge - fee"),
        plotter->get1DHistogram("cluster energy - bottom - no edge - fee"));

    for (int crystal_index = 2; crystal_index <= 6; ++crystal_index) { 
        
        ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - fee - crystal index = " + std::to_string(crystal_index)))->Divide(
            plotter->get1DHistogram("cluster count - matched - bottom - fee - crystal index = " + std::to_string(crystal_index)),
            plotter->get1DHistogram("cluster count - bottom - fee - crystal index = " + std::to_string(crystal_index)));
        
        ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - fee - crystal index = " + std::to_string(crystal_index)))->Divide(
            plotter->get1DHistogram("cluster count - matched - top - fee - crystal index = " + std::to_string(crystal_index)),
            plotter->get1DHistogram("cluster count - top - fee - crystal index = " + std::to_string(crystal_index)));
    }


    // Cluster time
    ((TGraphAsymmErrors*) plotter->buildGraph("cluster time - matched"))->Divide(
        plotter->get1DHistogram("cluster time - matched"),
        plotter->get1DHistogram("cluster time"));

    // Track plots
    /*
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - doca - top"))->Divide(
        plotter->get1DHistogram("doca - top - pass cuts"),
        plotter->get1DHistogram("doca - top"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - doca - bottom"))->Divide(
        plotter->get1DHistogram("doca - bottom - pass cuts"),
        plotter->get1DHistogram("doca - bottom"));

    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - z0 - top"))->Divide(
        plotter->get1DHistogram("z0 - top - pass cuts"),
        plotter->get1DHistogram("z0 - top"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - z0 - bottom"))->Divide(
        plotter->get1DHistogram("z0 - bottom - pass cuts"),
        plotter->get1DHistogram("z0 - bottom"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - sin(phi0) - top"))->Divide(
        plotter->get1DHistogram("sin(phi0) - top - pass cuts"),
        plotter->get1DHistogram("sin(phi0) - top"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - sin(phi0) - bottom"))->Divide(
        plotter->get1DHistogram("sin(phi0) - bottom - pass cuts"),
        plotter->get1DHistogram("sin(phi0) - bottom"));

    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - curvature - top"))->Divide(
        plotter->get1DHistogram("curvature - top - pass cuts"),
        plotter->get1DHistogram("curvature - top"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - curvature - bottom"))->Divide(
        plotter->get1DHistogram("curvature - bottom - pass cuts"),
        plotter->get1DHistogram("curvature - bottom"));

    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - tan_lambda - top"))->Divide(
        plotter->get1DHistogram("tan_lambda - top - pass cuts"),
        plotter->get1DHistogram("tan_lambda - top"));
    
    ((TGraphAsymmErrors*) plotter->buildGraph("tracking efficiency - tan_lambda - bottom"))->Divide(
        plotter->get1DHistogram("tan_lambda - bottom - pass cuts"),
        plotter->get1DHistogram("tan_lambda - bottom"));
    */

    plotter->saveToPdf("simple_tracking_efficiency.pdf");
    plotter->saveToRootFile("simple_tracking_efficiency.root");
}

void TrackClusterMatchingEfficiencyAnalysis::bookHistograms() { 

    // Plots for all clusters //
    ////////////////////////////

    // No cuts
    plotter->build1DHistogram("cluster energy", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");
    
    plotter->build2DHistogram("cluster energy v cluster size", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size")->GetYaxis()->SetTitle("Cluster size");
    
    plotter->build2DHistogram("cluster energy v cluster time", 50, 0, 1.5, 160, 0, 80);
    plotter->get2DHistogram("cluster energy v cluster time")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster time")->GetYaxis()->SetTitle("Cluster time (ns)");
    
    plotter->build2DHistogram("cluster count", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count")->GetYaxis()->SetTitle("Crystal Index - y");

    plotter->build2DHistogram("cluster position", 100, -200, 200, 50, -100, 100);
    plotter->get2DHistogram("cluster position")->GetXaxis()->SetTitle("Cluster Position - x (mm)");
    plotter->get2DHistogram("cluster position")->GetYaxis()->SetTitle("Cluster Position - y (mm)");

    plotter->build2DHistogram("cluster energy v cluster y", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y", 50, 0, 1.5, 12, -6, 6);
    
    // Plots of clusters split between top and bottom with no cuts
    plotter->build1DHistogram("cluster energy - top", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - top - fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - top", 160, 0, 80);
  
    plotter->build1DHistogram("cluster energy - bottom", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - bottom - fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - bottom", 160, 0, 80);

    // Time cut
    plotter->build2DHistogram("cluster energy v cluster time - time", 50, 0, 1.5, 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster y - time", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y - time", 50, 0, 1.5, 12, -6, 6);

    plotter->build2DHistogram("cluster energy v cluster size - time", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - time")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - time")->GetYaxis()->SetTitle("Cluster size");

    // Time + cluster size cut
    plotter->build2DHistogram("cluster energy v cluster time - time, size", 50, 0, 1.5, 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster y - time, size", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y - time, size", 50, 0, 1.5, 12, -6, 6);

    plotter->build2DHistogram("cluster energy v cluster size - time, size", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - time, size")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - time, size")->GetYaxis()->SetTitle("Cluster size");

    // Time + cluster size + no edge
    plotter->build2DHistogram("cluster energy v cluster time - time, size, no edge", 50, 0, 1.5, 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster y - time, size, no edge", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y - time, size, no edge", 50, 0, 1.5, 12, -6, 6);

    plotter->build2DHistogram("cluster energy v cluster size - time, size, no edge", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - time, size, no edge")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - time, size, no edge")->GetYaxis()->SetTitle("Cluster size");

    plotter->build1DHistogram("cluster energy - no edge", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - no edge", 160, 0, 80);

    // Plots of tracks //
    /////////////////////
    plotter->build1DHistogram("track time", 40, -20, 20)->GetXaxis()->SetTitle("Track time (ns)");
    plotter->build1DHistogram("E/p", 40, 0, 1.5)->GetXaxis()->SetTitle("E/p");
    plotter->build1DHistogram("E/p - no edge", 40, 0, 1.5)->GetXaxis()->SetTitle("E/p");

    // Track cluster matching plots //
    //////////////////////////////////

    plotter->build1DHistogram("cluster time - track time", 100, -20, 80);
    plotter->build2DHistogram("cluster x v extrapolated track x", 100, -100, 100, 100, -100, 100);
    plotter->build2DHistogram("cluster y v extrapolated track y", 100, -100, 100, 100, -100, 100);
    plotter->build1DHistogram("cluster x - extrapolated track x", 50, -25, 25);
    plotter->build1DHistogram("cluster y - extrapolated track y", 50, -25, 25);

    // Plots of tracks matched to clusters

    plotter->build2DHistogram("cluster count - matched", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - matched")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - matched")->GetYaxis()->SetTitle("Crystal Index - y");

    plotter->build1DHistogram("cluster energy - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");
    plotter->build1DHistogram("E/p - matched", 40, 0, 1.5)->GetXaxis()->SetTitle("E/p");

    // Plots of tracks matched to clusters split into top and bottom tracks
    plotter->build1DHistogram("cluster energy - matched - top", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched - top", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster energy - matched - bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched - bottom", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    // Plots of FEE //
    //////////////////

    plotter->build2DHistogram("cluster count - fee", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - fee")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - fee")->GetYaxis()->SetTitle("Crystal Index - y");
    
    plotter->build2DHistogram("cluster count - matched - fee", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - matched - fee")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - matched - fee")->GetYaxis()->SetTitle("Crystal Index - y");
    
    plotter->build1DHistogram("E/p - matched - fee", 40, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - fee", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
 
    plotter->build2DHistogram("cluster energy v cluster time - fee", 50, 0, 1.5, 160, 0, 80);
    plotter->get2DHistogram("cluster energy v cluster time - fee")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster time - fee")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster energy v cluster size - fee", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - fee")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - fee")->GetYaxis()->SetTitle("Cluster size");

    // FEE plots split into top and bottom
    plotter->build1DHistogram("cluster energy - matched - top - fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - bottom - fee", 50, 0, 1.5);
    
    plotter->build1DHistogram("cluster energy - top - no edge - fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - bottom - no edge - fee", 50, 0, 1.5);
    
    for (int crystal_index = 2; crystal_index <= 6; ++crystal_index) { 
        plotter->build1DHistogram("cluster count - top - fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");
        
        plotter->build1DHistogram("cluster count - bottom - fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");

        plotter->build1DHistogram("cluster count - matched - top - fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");
        
        plotter->build1DHistogram("cluster count - matched - bottom - fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");

    }

    plotter->build1DHistogram("cluster energy - matched - top - no edge - fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - bottom - no edge - fee", 50, 0, 1.5);

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

    plotter->build1DHistogram("doca - top", 80, -10, 10);
    plotter->build1DHistogram("z0 - top", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - top", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - top", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - top", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - top", 40, -0.1, 0.1);

    plotter->build1DHistogram("doca - top - pass cuts", 80, -10, 10);
    plotter->build1DHistogram("z0 - top - pass cuts", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - top - pass cuts", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - top - pass cuts", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - top - pass cuts", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - top - pass cuts", 40, -0.1, 0.1);

    plotter->build1DHistogram("doca - bottom", 80, -10, 10);
    plotter->build1DHistogram("z0 - bottom", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - bottom", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - bottom", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - bottom", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - bottom", 40, -0.1, 0.1);

    plotter->build1DHistogram("doca - bottom - pass cuts", 80, -10, 10);
    plotter->build1DHistogram("z0 - bottom - pass cuts", 80, -2, 2);
    plotter->build1DHistogram("sin(phi0) - bottom - pass cuts", 40, -0.2, 0.2);
    plotter->build1DHistogram("curvature - bottom - pass cuts", 50, -0.001, 0.001);
    plotter->build1DHistogram("tan_lambda - bottom - pass cuts", 100, -0.1, 0.1);
    plotter->build1DHistogram("cos(theta) - bottom - pass cuts", 40, -0.1, 0.1);
    
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

std::string TrackClusterMatchingEfficiencyAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;
}


bool TrackClusterMatchingEfficiencyAnalysis::passEnergyCut(EcalCluster* cluster) { 
    if (cluster->getEnergy() < cluster_energy_low_threshold
            || cluster->getEnergy() > cluster_energy_high_threshold) return false;

    return true;
}

bool TrackClusterMatchingEfficiencyAnalysis::passClusterTimeCut(EcalCluster* cluster) {   
    if (cluster->getClusterTime() < 41 || cluster->getClusterTime() > 50) return false;

    return true;   
}

bool TrackClusterMatchingEfficiencyAnalysis::passClusterSizeCut(EcalCluster* cluster) { 
    if (cluster->getEcalHits()->GetEntriesFast() < 3) return false;
    return true;
}

bool TrackClusterMatchingEfficiencyAnalysis::isEdgeCrystal(EcalHit* hit) { 
    
    int y_crystal_index = hit->getYCrystalIndex();

    if (abs(y_crystal_index) == 1) return true;
    return false;
}


bool TrackClusterMatchingEfficiencyAnalysis::isMatch(EcalCluster* cluster, SvtTrack* track) { 

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

