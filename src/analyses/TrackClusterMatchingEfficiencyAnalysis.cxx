
#include <TrackClusterMatchingEfficiencyAnalysis.h>

TrackClusterMatchingEfficiencyAnalysis::TrackClusterMatchingEfficiencyAnalysis() 
    : track(NULL),
      cluster(NULL),
      plotter(new Plotter()),
      matcher(new TrackClusterMatcher()),
      cluster_energy_low_threshold(.7 /* GeV */),
      cluster_energy_high_threshold(1.15 /* GeV */),
      cuts_enabled(true),
      class_name("TrackClusterMatchingEfficiencyAnalysis"),
      total_events(0),
      total_single1_triggers(0),
      total_events_with_tracks(0) {

}

TrackClusterMatchingEfficiencyAnalysis::~TrackClusterMatchingEfficiencyAnalysis() {
    delete plotter;
    delete matcher;
}

void TrackClusterMatchingEfficiencyAnalysis::initialize() { 
    this->bookHistograms();
    
    matcher->enablePlots(true);
}

void TrackClusterMatchingEfficiencyAnalysis::processEvent(HpsEvent* event) { 

    // Increment the total events counter
    total_events++; 

    // Only look at single 1 triggers
    if (!event->isSingle1Trigger()) return;
    
    // Increment the singles1 trigger counter
    total_single1_triggers++;

    if (event->getNumberOfTracks() != 0) total_events_with_tracks++;
    
    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    matched_tracks.clear();
    // Loop over all clusters and apply cuts to try and isolate FEE's using
    // just Ecal information 
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
  
        // Get the cluster from the event 
        cluster = event->getEcalCluster(cluster_n);
  
        // Get the energy associated with this cluster 
        double cluster_energy = cluster->getEnergy();

        //  Get the time associated with the cluster
        double cluster_time = cluster->getClusterTime();

        // Fill the cluster information for all events
        plotter->get1DHistogram("cluster energy")->Fill(cluster_energy);
        plotter->get1DHistogram("cluster time")->Fill(cluster_time);
        plotter->get2DHistogram("cluster energy v cluster time")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster position")->Fill(
                cluster->getPosition()[0], cluster->getPosition()[1]);

        // Get the seed hit of the cluster
        EcalHit* seed_hit = cluster->getSeed();

        plotter->get2DHistogram("cluster energy v cluster seed energy")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        plotter->get2DHistogram("cluster count")->Fill(
                seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);

        // Make the same plots for the top and bottom Ecal volumes. Use the
        // seed crystal index to distinguish between cluster that are in the top
        // or bottom Ecal volume 
        if (seed_hit->getYCrystalIndex() > 0) { 
            plotter->get1DHistogram("cluster energy - top")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - top")->Fill(cluster_time);
        } else { 
            plotter->get1DHistogram("cluster energy - bottom")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - bottom")->Fill(cluster_time);
        }  
    
        // Before applying any cuts, check if the cluster is matched to a track.  This will
        // allow calculation of the track-cluster matching efficiencies as a function
        // of energy and time.
        if (matcher->getMatchingTrack(cluster) != NULL) { 
        
            plotter->get2DHistogram("cluster count - matched")->Fill(
                    seed_hit->getXCrystalIndex(), 
                    seed_hit->getYCrystalIndex(), 1);
       
            plotter->get1DHistogram("cluster energy - matched")->Fill(cluster->getEnergy(), 1);
            plotter->get1DHistogram("cluster time - matched")->Fill(cluster->getClusterTime(), 1);

            if (seed_hit->getYCrystalIndex() > 0) { 
                plotter->get1DHistogram("cluster energy - matched - top")->Fill(cluster->getEnergy(), 1);
                plotter->get1DHistogram("cluster time - matched - top")->Fill(cluster->getClusterTime(), 1);
            } else { 
                plotter->get1DHistogram("cluster energy - matched - bottom")->Fill(cluster->getEnergy(), 1);
                plotter->get1DHistogram("cluster time - matched - bottom")->Fill(cluster->getClusterTime(), 1);
            }  
        }

        if (!isEdgeCrystal(seed_hit)) { 
            plotter->get1DHistogram("cluster energy - no edge")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - no edge")->Fill(cluster_time);
        }
 
        // Check that the cluster passes the time requirement
        if (!passClusterTimeCut(cluster)) return;

        plotter->get1DHistogram("cluster time - cuts: time")->Fill(cluster_time);
        plotter->get2DHistogram("cluster energy v cluster time - cuts: time")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster y - cuts: time")->Fill(cluster_energy, 
                cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster energy v crystal index - y - cuts: time")->Fill(cluster_energy,
                seed_hit->getYCrystalIndex());
        plotter->get2DHistogram("cluster energy v cluster size - cuts: time")->Fill(cluster_energy, 
                    cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        // Check that the cluster passes the size requirement
        if (!passClusterSizeCut(cluster)) return; 

        plotter->get2DHistogram("cluster energy v cluster time - cuts: time, size")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster y - cuts: time, size")->Fill(cluster_energy, cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster energy v crystal index - y - cuts: time, size")->Fill(cluster_energy, seed_hit->getYCrystalIndex());
        plotter->get2DHistogram("cluster energy v cluster size - cuts: time, size")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        // Check that the cluster passes the energy requirement
        if (!passEnergyCut(cluster)) return;
    
        plotter->get2DHistogram("cluster energy v cluster time - cuts: fee")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size - cuts: fee")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster count - cuts: fee")->Fill(
                seed_hit->getXCrystalIndex(), seed_hit->getYCrystalIndex(), 1);
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: fee")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        int crystal_index = (int) abs(seed_hit->getYCrystalIndex());
        if (seed_hit->getYCrystalIndex() > 0) { 
            plotter->get1DHistogram("cluster energy - top - cuts: fee")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster count - top - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
            if (!isEdgeCrystal(seed_hit)) {     
                plotter->get1DHistogram("cluster energy - top - no edge - cuts: fee")->Fill(cluster_energy);
            }

        } else { 
            plotter->get1DHistogram("cluster energy - bottom - cuts: fee")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster count - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
            if (!isEdgeCrystal(seed_hit)) { 
                plotter->get1DHistogram("cluster energy - bottom - no edge - cuts: fee")->Fill(cluster_energy);
            }
        } 

        // Check if the cluster has a matching track
        if (matcher->getMatchingTrack(cluster) != NULL) { 
            
            track = matcher->getMatchingTrack(cluster); 

            // Calculate the momentum magnitude 
            std::vector<double> p = track->getMomentum();
            double p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
            double pt = sqrt(p[0]*p[0] + p[1]*p[1]);

            plotter->get2DHistogram("cluster count - matched - cuts: fee")->Fill(
                    seed_hit->getXCrystalIndex(), 
                    seed_hit->getYCrystalIndex(), 1);

            int crystal_index = (int) abs(seed_hit->getYCrystalIndex());
            if (seed_hit->getYCrystalIndex() > 0) { 
                plotter->get1DHistogram("cluster energy - matched - top - cuts: fee")->Fill(cluster_energy);
                plotter->get1DHistogram("cluster count - matched - top - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                if (!isEdgeCrystal(seed_hit)) {   
                    plotter->get1DHistogram("cluster energy - matched - top - cuts: fee, no edge")->Fill(cluster_energy);
                    plotter->get1DHistogram("cluster count - matched - top - no edge - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                }  
            } else { 
                plotter->get1DHistogram("cluster energy - matched - bottom - cuts: fee")->Fill(cluster_energy);
                plotter->get1DHistogram("cluster count - matched - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                if (!isEdgeCrystal(seed_hit)) {   
                    plotter->get1DHistogram("cluster energy - matched - bottom - cuts: fee")->Fill(cluster_energy);
                    plotter->get1DHistogram("cluster count - matched - bottom - no edge - cuts: fee - crystal index = " + std::to_string(crystal_index))->Fill(seed_hit->getXCrystalIndex(), 1);
                }  
            }
        
            if (track->isTopTrack()) { 
                plotter->get1DHistogram("p - matched - top - cuts: fee")->Fill(p_mag);
                plotter->get1DHistogram("pt - matched - top - cuts: fee")->Fill(pt);
                plotter->get1DHistogram("px - matched - top - cuts: fee")->Fill(p[0]);
                plotter->get1DHistogram("py - matched - top - cuts: fee")->Fill(p[1]);
                plotter->get1DHistogram("pz - matched - top - cuts: fee")->Fill(p[2]);
            } else if (track->isBottomTrack()) { 
                plotter->get1DHistogram("p - matched - bottom - cuts: fee")->Fill(p_mag);
                plotter->get1DHistogram("pt - matched - bottom - cuts: fee")->Fill(pt);
                plotter->get1DHistogram("px - matched - bottom - cuts: fee")->Fill(p[0]);
                plotter->get1DHistogram("py - matched - bottom - cuts: fee")->Fill(p[1]);
                plotter->get1DHistogram("pz - matched - bottom - cuts: fee")->Fill(p[2]);
            }
        
        } else { 
            
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
    plotter->get2DHistogram("cluster count - matched - cuts: fee")->Divide(
            plotter->get2DHistogram("cluster count - cuts: fee"));

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

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - cuts: fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - top - cuts: fee"),
        plotter->get1DHistogram("cluster energy - top - cuts: fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - cuts: fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - bottom - cuts: fee"),
        plotter->get1DHistogram("cluster energy - bottom - cuts: fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - no edge - cuts: fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - top - no edge - cuts: fee"),
        plotter->get1DHistogram("cluster energy - top - no edge - cuts: fee"));

    ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - no edge - cuts: fee"))->Divide(
        plotter->get1DHistogram("cluster energy - matched - bottom - no edge - cuts: fee"),
        plotter->get1DHistogram("cluster energy - bottom - no edge - cuts: fee"));

    for (int crystal_index = 2; crystal_index <= 6; ++crystal_index) { 
        
        ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index)))->Divide(
            plotter->get1DHistogram("cluster count - matched - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index)),
            plotter->get1DHistogram("cluster count - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index)));
        
        ((TGraphAsymmErrors*) plotter->buildGraph("cluster-track match efficiency - top - cuts: fee - crystal index = " + std::to_string(crystal_index)))->Divide(
            plotter->get1DHistogram("cluster count - matched - top - cuts: fee - crystal index = " + std::to_string(crystal_index)),
            plotter->get1DHistogram("cluster count - top - cuts: fee - crystal index = " + std::to_string(crystal_index)));
    }


    // Cluster time
    ((TGraphAsymmErrors*) plotter->buildGraph("cluster time - matched"))->Divide(
        plotter->get1DHistogram("cluster time - matched"),
        plotter->get1DHistogram("cluster time"));

    plotter->saveToPdf("simple_tracking_efficiency.pdf");
    plotter->saveToRootFile("simple_tracking_efficiency.root");

    matcher->saveHistograms();
}

void TrackClusterMatchingEfficiencyAnalysis::bookHistograms() { 

    // Plots for all clusters //
    ////////////////////////////

    // No cuts
    plotter->build1DHistogram("cluster energy", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");
   
    plotter->build2DHistogram("cluster energy v cluster seed energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster energy v cluster seed energy")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster seed energy")->GetYaxis()->SetTitle("Cluster seed energy (GeV)");
    
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
    plotter->build1DHistogram("cluster energy - top - cuts: fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - top", 160, 0, 80);
  
    plotter->build1DHistogram("cluster energy - bottom", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - bottom - cuts: fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - bottom", 160, 0, 80);

    // Time cut
    plotter->build1DHistogram("cluster time - cuts: time", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");
    
    plotter->build2DHistogram("cluster energy v cluster seed energy - cuts: time", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time")->GetYaxis()->SetTitle("Cluster seed energy (GeV)");

    plotter->build2DHistogram("cluster energy v cluster time - cuts: time", 50, 0, 1.5, 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster y - cuts: time", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y - cuts: time", 50, 0, 1.5, 12, -6, 6);

    plotter->build2DHistogram("cluster energy v cluster size - cuts: time", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - cuts: time")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - cuts: time")->GetYaxis()->SetTitle("Cluster size");

    // Time + cluster size cut
    plotter->build2DHistogram("cluster energy v cluster time - cuts: time, size", 50, 0, 1.5, 160, 0, 80);

    plotter->build2DHistogram("cluster energy v cluster y - cuts: time, size", 50, 0, 1.5, 50, -100, 100);

    plotter->build2DHistogram("cluster energy v crystal index - y - cuts: time, size", 50, 0, 1.5, 12, -6, 6);

    plotter->build2DHistogram("cluster energy v cluster seed energy - cuts: time, size", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size")->GetYaxis()->SetTitle("Cluster seed energy (GeV)");

    plotter->build2DHistogram("cluster energy v cluster size - cuts: time, size", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - cuts: time, size")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - cuts: time, size")->GetYaxis()->SetTitle("Cluster size");

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

    // Plots of tracks matched to clusters

    plotter->build2DHistogram("cluster count - matched", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - matched")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - matched")->GetYaxis()->SetTitle("Crystal Index - y");

    plotter->build1DHistogram("cluster energy - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    // Plots of tracks matched to clusters split into top and bottom tracks
    plotter->build1DHistogram("cluster energy - matched - top", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched - top", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    plotter->build1DHistogram("cluster energy - matched - bottom", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched - bottom", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    // Plots of FEE //
    //////////////////

    plotter->build2DHistogram("cluster energy v cluster time - cuts: fee", 50, 0, 1.5, 160, 0, 80);
    plotter->get2DHistogram("cluster energy v cluster time - cuts: fee")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster time - cuts: fee")->GetYaxis()->SetTitle("Cluster time (ns)");

    plotter->build2DHistogram("cluster energy v cluster size - cuts: fee", 50, 0, 1.5, 10, 0, 10);
    plotter->get2DHistogram("cluster energy v cluster size - cuts: fee")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster size - cuts: fee")->GetYaxis()->SetTitle("Cluster size");

    plotter->build2DHistogram("cluster count - cuts: fee", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - cuts: fee")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - cuts: fee")->GetYaxis()->SetTitle("Crystal Index - y");
   
    plotter->build2DHistogram("cluster energy v cluster seed energy - cuts: fee", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: fee")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: fee")->GetYaxis()->SetTitle("Cluster seed energy (GeV)");

    // 
    plotter->build2DHistogram("cluster count - matched - cuts: fee", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - matched - cuts: fee")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - matched - cuts: fee")->GetYaxis()->SetTitle("Crystal Index - y");
    
    plotter->build1DHistogram("cluster energy - matched - cuts: fee", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");

    // FEE plots split into top and bottom
    plotter->build1DHistogram("cluster energy - matched - top - cuts: fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - bottom - cuts: fee", 50, 0, 1.5);
    
    plotter->build1DHistogram("cluster energy - matched - top - no edge - cuts: fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - bottom - no edge - cuts: fee", 50, 0, 1.5);

    plotter->build1DHistogram("cluster energy - top - no edge - cuts: fee", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - bottom - no edge - cuts: fee", 50, 0, 1.5);
    
    for (int crystal_index = 1; crystal_index <= 6; ++crystal_index) { 
        plotter->build1DHistogram("cluster count - top - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");
        
        plotter->build1DHistogram("cluster count - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");

        plotter->build1DHistogram("cluster count - matched - top - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");
        
        plotter->build1DHistogram("cluster count - matched - bottom - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");


        plotter->build1DHistogram("cluster count - matched - top - no edge - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");
        
        plotter->build1DHistogram("cluster count - matched - bottom - no edge - cuts: fee - crystal index = " + std::to_string(crystal_index),
                47, -23, 24)->GetXaxis()->SetTitle("Crystal Index - x");

    }

    plotter->build1DHistogram("cluster energy - matched - top - cuts: fee, no edge", 50, 0, 1.5);
    plotter->build1DHistogram("cluster energy - matched - bottom - cuts: fee, no edge", 50, 0, 1.5);

    plotter->build1DHistogram("p - matched - top - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - matched - top - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - matched - top - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - matched - top - cuts: fee", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - matched - top - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    plotter->build1DHistogram("p - matched - bottom - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - matched - bottom - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - matched - bottom - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - matched - bottom - cuts: fee", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - matched - bottom - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    plotter->build2DHistogram("cluster count - no match", 47, -23, 24, 12, -6, 6);

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
    if (cluster->getClusterTime() < 41.5 || cluster->getClusterTime() > 49.67) return false;

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

 
    // Check that dx and dy between the extrapolated track and cluster
    // positions is reasonable
    if ((track->isTopTrack() && (delta_x > 12. || delta_x < -7.61)) ||
        (track->isBottomTrack() && (delta_x > 6.0 || delta_x < -13.75))) return false;

    if ((track->isTopTrack() && (delta_y > 14 || delta_y < -14)) ||
        (track->isBottomTrack() && (delta_y > 14 || delta_y < -14))) return false;

    //if (cluster->getEnergy()/p < .5) return false;

    return true;
}
