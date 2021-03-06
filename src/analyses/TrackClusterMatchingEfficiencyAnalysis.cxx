
#include <TrackClusterMatchingEfficiencyAnalysis.h>

TrackClusterMatchingEfficiencyAnalysis::TrackClusterMatchingEfficiencyAnalysis() 
    : cluster(NULL),
      plotter(new Plotter()),
      matcher(new TrackClusterMatcher()),
      cluster_energy_low_threshold(.8 /* GeV */),
      cluster_energy_high_threshold(1.1 /* GeV */),
      cuts_enabled(true),
      class_name("TrackClusterMatchingEfficiencyAnalysis"),
      event_counter(0),
      event_track_counter(0), 
      bias_on_counter(0), 
      single1_trigger_counter(0),
      ecal_cluster_time_cut_pass_counter(0), 
      ecal_cluster_size_cut_pass_counter(0), 
      ecal_cluster_energy_cut_pass_counter(0), 
      ecal_cluster_seed_energy_cut_pass_counter(0), 
      svt_closed_position_counter(0) {
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
    event_counter++; 

    // Only look at single 1 triggers
    if (false) { 
        if (!event->isSingle1Trigger()) return;


        // Only look at events with the SVT bias ON
        if (!event->isSvtBiasOn()) return; 
    
        // Increment the counter keeping track of events with SVT bias ON
        bias_on_counter++; 

        // Only look at events where the SVT is closed
        if (!event->isSvtClosed()) return;

        // Increment the counter keeping track of events with the SVT bias ON and
        // the SVT closed 
        svt_closed_position_counter++; 
    }

        // Increment the singles1 trigger counter
        single1_trigger_counter++;

    if (event->getNumberOfTracks() != 0) event_track_counter++;
    
    // Find all track-cluster matches in the event
    matcher->findAllMatches(event);

    // Loop over all clusters and apply cuts to try and isolate FEE's using
    // just Ecal information 
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) {
  
        // Get the cluster from the event 
        cluster = event->getEcalCluster(cluster_n);
  
        // Get the energy associated with this cluster 
        double cluster_energy = cluster->getEnergy();

        //  Get the time associated with the cluster
        double cluster_time = cluster->getClusterTime();

        // Get the size of the cluster
        int cluster_size = cluster->getEcalHits()->GetEntriesFast(); 

        // Fill the cluster information for all events
        plotter->get1DHistogram("cluster energy")->Fill(cluster_energy);
        plotter->get1DHistogram("cluster time")->Fill(cluster_time);
        if (cluster_size < 10) { 
            plotter->get1DHistogram("cluster time - cluster size " + std::to_string(cluster_size))->Fill(cluster_time);
        }
        plotter->get2DHistogram("cluster energy v cluster time")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size")->Fill(cluster_energy, cluster_size);
        plotter->get2DHistogram("cluster position")->Fill(cluster->getPosition()[0], cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster time v cluster size")->Fill(cluster_time, cluster_size); 

        // Get the seed hit of the cluster
        EcalHit* seed_hit = cluster->getSeed();
        
        // Get the crystal indices of the cluster
        int crystal_index_x = seed_hit->getXCrystalIndex(); 
        int crystal_index_y = seed_hit->getYCrystalIndex(); 

        plotter->get2DHistogram("cluster energy v cluster seed energy")->Fill(cluster_energy, seed_hit->getEnergy());
        plotter->get2DHistogram("cluster count")->Fill(crystal_index_x, crystal_index_y, 1);

        // Make the same plots for the top and bottom Ecal volumes. Use the
        // seed crystal index to distinguish between cluster that are in the top
        // or bottom Ecal volume 
        if (crystal_index_y > 0) { 
            plotter->get1DHistogram("cluster energy - top")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - top")->Fill(cluster_time);
        } else { 
            plotter->get1DHistogram("cluster energy - bottom")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - bottom")->Fill(cluster_time);
        }  
    
        // Get the track was matched to this cluster
        SvtTrack* track = matcher->getMatchingTrack(cluster); 

        // Before applying any cuts, check if the cluster is matched to a track.  This will
        // allow calculation of the track-cluster matching efficiencies as a function
        // of energy and time. 
        if (track != NULL) { 
        
            plotter->get2DHistogram("cluster count - matched")->Fill(crystal_index_x, crystal_index_y, 1);
            plotter->get1DHistogram("cluster energy - matched")->Fill(cluster_energy, 1);
            plotter->get1DHistogram("cluster time - matched")->Fill(cluster_time, 1);

            if (crystal_index_y > 0) { 
                plotter->get1DHistogram("cluster energy - matched - top")->Fill(cluster_energy, 1);
                plotter->get1DHistogram("cluster time - matched - top")->Fill(cluster_time, 1);
            } else { 
                plotter->get1DHistogram("cluster energy - matched - bottom")->Fill(cluster_energy, 1);
                plotter->get1DHistogram("cluster time - matched - bottom")->Fill(cluster_time, 1);
            }  
        }

        if (!isEdgeCrystal(seed_hit)) { 
            plotter->get1DHistogram("cluster energy - no edge")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster time - no edge")->Fill(cluster_time);
        }
 
        // Check that the cluster passes the time requirement
        if (!passClusterTimeCut(cluster)) return;
        ecal_cluster_time_cut_pass_counter++; 

        plotter->get1DHistogram("cluster time - cuts: time")->Fill(cluster_time);
        plotter->get2DHistogram("cluster energy v cluster time - cuts: time")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster y - cuts: time")->Fill(cluster_energy, 
                cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster energy v crystal index - y - cuts: time")->Fill(cluster_energy,
                crystal_index_y);
        plotter->get2DHistogram("cluster energy v cluster size - cuts: time")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        // Check that the cluster passes the size requirement
        if (!passClusterSizeCut(cluster)) return; 
        ecal_cluster_size_cut_pass_counter++; 

        plotter->get2DHistogram("cluster energy v cluster time - cuts: time, size")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster y - cuts: time, size")->Fill(cluster_energy, 
                cluster->getPosition()[1]);
        plotter->get2DHistogram("cluster energy v crystal index - y - cuts: time, size")->Fill(cluster_energy, 
                crystal_index_y);
        plotter->get2DHistogram("cluster energy v cluster size - cuts: time, size")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        if (!isEdgeCrystal(seed_hit)) { 
            plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size, no edge")->Fill(
                    cluster_energy, seed_hit->getEnergy());
        }

        // Check that the cluster passes the energy requirement
        if (!passEnergyCut(cluster)) return;
        ecal_cluster_energy_cut_pass_counter++; 

        if (seed_hit->getEnergy() < .4) return; 
        ecal_cluster_seed_energy_cut_pass_counter++; 

        plotter->get2DHistogram("cluster energy v cluster time - cuts: fee")->Fill(cluster_energy, cluster_time);
        plotter->get2DHistogram("cluster energy v cluster size - cuts: fee")->Fill(cluster_energy, 
                cluster->getEcalHits()->GetEntriesFast());
        plotter->get2DHistogram("cluster count - cuts: fee")->Fill(crystal_index_x, crystal_index_y, 1);
        plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: fee")->Fill(cluster_energy, 
                seed_hit->getEnergy());

        int crystal_index = (int) abs(crystal_index_y);
        if (crystal_index_y > 0) { 
            plotter->get1DHistogram("cluster energy - top - cuts: fee")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster count - top - cuts: fee - crystal index = " 
                    + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
            if (!isEdgeCrystal(seed_hit)) {     
                plotter->get1DHistogram("cluster energy - top - no edge - cuts: fee")->Fill(cluster_energy);
            }

        } else { 
            plotter->get1DHistogram("cluster energy - bottom - cuts: fee")->Fill(cluster_energy);
            plotter->get1DHistogram("cluster count - bottom - cuts: fee - crystal index = "
                    + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
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

            plotter->get2DHistogram("cluster count - matched - cuts: fee")->Fill(crystal_index_x, crystal_index_y, 1);
            plotter->get1DHistogram("cluster energy - matched - cuts: fee")->Fill(cluster_energy);
            plotter->get2DHistogram("cluster energy v p - matched - cuts: fee")->Fill(p_mag, cluster_energy); 
            
            int crystal_index = (int) abs(crystal_index_y);
            if (crystal_index_y > 0) { 
                plotter->get1DHistogram("cluster energy - matched - top - cuts: fee")->Fill(cluster_energy);
                plotter->get1DHistogram("cluster count - matched - top - cuts: fee - crystal index = " 
                        + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
                if (!isEdgeCrystal(seed_hit)) {   
                    plotter->get1DHistogram("cluster energy - matched - top - no edge - cuts: fee")->Fill(
                            cluster_energy);
                    plotter->get1DHistogram("cluster count - matched - top - no edge - cuts: fee - crystal index = " 
                            + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
                }  
            } else { 
                plotter->get1DHistogram("cluster energy - matched - bottom - cuts: fee")->Fill(cluster_energy);
                plotter->get1DHistogram("cluster count - matched - bottom - cuts: fee - crystal index = " 
                        + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
                if (!isEdgeCrystal(seed_hit)) {   
                    plotter->get1DHistogram("cluster energy - matched - bottom - no edge - cuts: fee")->Fill(
                            cluster_energy);
                    plotter->get1DHistogram("cluster count - matched - bottom - no edge - cuts: fee - crystal index = " 
                            + std::to_string(crystal_index))->Fill(crystal_index_x, 1);
                }  
            }
      
            if (!isEdgeCrystal(seed_hit)) { 

            // Loop over all of the tracks in the event
            GblTrack* gbl_track = NULL; 
            for (int gbl_track_n = 0; gbl_track_n < event->getNumberOfGblTracks(); ++gbl_track_n) { 

                // Get a GBL track from the event
                //gbl_track = event->getGblTrack(gbl_track_n); 
        
                // Get the seed track associated with the GBL track
                //SvtTrack* seed_track = (SvtTrack*) gbl_track->getSeedTrack().GetObject();
            
                //if (seed_track == track) break;  
            }
            
            std::vector<double> gbl_p;
            double gbl_p_mag;
            double gbl_pt;
            if (gbl_track != NULL) { 
                gbl_p = gbl_track->getMomentum();
                gbl_p_mag = sqrt(gbl_p[0]*gbl_p[0] + gbl_p[1]*gbl_p[1] + gbl_p[2]*gbl_p[2]);
                gbl_pt = sqrt(gbl_p[0]*gbl_p[0] + gbl_p[1]*gbl_p[1]);
            }

            if (track->isTopTrack()) { 
                plotter->get1DHistogram("p - matched - top - cuts: fee")->Fill(p_mag);
                plotter->get1DHistogram("pt - matched - top - cuts: fee")->Fill(pt);
                plotter->get1DHistogram("px - matched - top - cuts: fee")->Fill(p[0]);
                plotter->get1DHistogram("py - matched - top - cuts: fee")->Fill(p[1]);
                plotter->get1DHistogram("pz - matched - top - cuts: fee")->Fill(p[2]);

                plotter->get1DHistogram("doca - matched - top")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - matched - top")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - matched - top")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - matched - top")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - matched - top")->Fill(track->getTanLambda()); 
                plotter->get1DHistogram("cos(theta) - matched - top")->Fill(TrackExtrapolator::getCosTheta(track));

                if (gbl_track != NULL) { 
                    plotter->get1DHistogram("p - matched - gbl - top - cuts: fee")->Fill(gbl_p_mag);
                    plotter->get1DHistogram("doca - matched - gbl - top")->Fill(gbl_track->getD0());
                    plotter->get1DHistogram("z0 - matched - gbl - top")->Fill(gbl_track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - matched - gbl - top")->Fill(sin(gbl_track->getPhi0()));
                    plotter->get1DHistogram("curvature - matched - gbl - top")->Fill(gbl_track->getOmega());
                    //plotter->get1DHistogram("cos(theta) - matched - gbl - top")->Fill(TrackExtrapolator::getCosTheta(gbl_track));
                
                }

            } else if (track->isBottomTrack()) { 
                plotter->get1DHistogram("p - matched - bottom - cuts: fee")->Fill(p_mag);
                plotter->get1DHistogram("pt - matched - bottom - cuts: fee")->Fill(pt);
                plotter->get1DHistogram("px - matched - bottom - cuts: fee")->Fill(p[0]);
                plotter->get1DHistogram("py - matched - bottom - cuts: fee")->Fill(p[1]);
                plotter->get1DHistogram("pz - matched - bottom - cuts: fee")->Fill(p[2]);

                plotter->get1DHistogram("doca - matched - bottom")->Fill(track->getD0());
                plotter->get1DHistogram("z0 - matched - bottom")->Fill(track->getZ0());
                plotter->get1DHistogram("sin(phi0) - matched - bottom")->Fill(sin(track->getPhi0()));
                plotter->get1DHistogram("curvature - matched - bottom")->Fill(track->getOmega());
                plotter->get1DHistogram("tan_lambda - matched - bottom")->Fill(track->getTanLambda()); 
                plotter->get1DHistogram("cos(theta) - matched - bottom")->Fill(TrackExtrapolator::getCosTheta(track));

                if (gbl_track != NULL) { 
                    plotter->get1DHistogram("p - matched - gbl - bottom - cuts: fee")->Fill(gbl_p_mag);
                    plotter->get1DHistogram("doca - matched - gbl - bottom")->Fill(gbl_track->getD0());
                    plotter->get1DHistogram("z0 - matched - gbl - bottom")->Fill(gbl_track->getZ0());
                    plotter->get1DHistogram("sin(phi0) - matched - gbl - bottom")->Fill(sin(gbl_track->getPhi0()));
                    plotter->get1DHistogram("curvature - matched - gbl - bottom")->Fill(gbl_track->getOmega());
                    //plotter->get1DHistogram("cos(theta) - matched - gbl - bottom")->Fill(TrackExtrapolator::getCosTheta(gbl_track));
                }
            }
        }
        
        } else { 
            
            plotter->get2DHistogram("cluster count - no match")->Fill(
                    crystal_index_x, crystal_index_y, 1);
        }
    }
}

void TrackClusterMatchingEfficiencyAnalysis::finalize() { 

    std::cout << "//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//" << std::endl;  
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of events: " << event_counter << std::endl;
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of single1 triggers: " 
              << single1_trigger_counter << std::endl;
    /*std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of single1 trigger events with SVT "
              << "bias ON: " << bias_on_counter << std::endl; 
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of single1 trigger events with SVT "
              << "bias ON and in closed position: " << svt_closed_position_counter << std::endl; */
    std::cout << "[ TrackClusterMatchingEfficiencyAnalysis ] Total number of events with tracks: "
              << event_track_counter << std::endl;
    std::cout << "Pass cluster time cut: " << ecal_cluster_time_cut_pass_counter << std::endl;
    std::cout << "Pass cluster size cut: " << ecal_cluster_size_cut_pass_counter << std::endl;
    std::cout << "Pass cluster energy cut: " << ecal_cluster_energy_cut_pass_counter << std::endl;
    std::cout << "Pass cluster seed cut: " << ecal_cluster_seed_energy_cut_pass_counter << std::endl;
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

    for (int crystal_index = 2; crystal_index < 6; ++crystal_index) { 
        
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

    // Fit the momentum distributions
    TFile* roo_fits_file = new TFile("fits.root", "RECREATE"); 
    RooPlot* plot = NULL;

    RooRealVar p_var_top("p_var_top", "Top Track FEE Momentum (GeV)", .2, 2.0); 
    plot = RooFitter::fitToGaussian(plotter->get1DHistogram("p - matched - top - cuts: fee"), p_var_top);
    //plot = RooFitter::fitToDoubleGaussian(plotter->get1DHistogram("p - matched - top - cuts: fee"), p_var);
    plot->Draw(); 
    plot->Write();
    //plot->pullHist() 

    RooRealVar p_var_bot("p_var_bot", "Bottom Track FEE Momentum (GeV)", .2, 2.0); 
    plot = RooFitter::fitToGaussian(plotter->get1DHistogram("p - matched - bottom - cuts: fee"), p_var_bot);
    //plot = RooFitter::fitToDoubleGaussian(plotter->get1DHistogram("p - matched - bottom - cuts: fee"), p_var);
    plot->Draw(); 
    plot->Write(); 

    //plot = RooFitter::fitToGaussian(plotter->get1DHistogram("p - matched - gbl - top - cuts: fee"), p_var);
    //plot = RooFitter::fitToDoubleGaussian(plotter->get1DHistogram("p - matched - gbl - top - cuts: fee"), p_var);
    plot->Draw(); 
    plot->Write();
    //plot->pullHist() 

    //plot = RooFitter::fitToGaussian(plotter->get1DHistogram("p - matched - gbl - bottom - cuts: fee"), p_var);
    //plot = RooFitter::fitToDoubleGaussian(plotter->get1DHistogram("p - matched - gbl - bottom - cuts: fee"), p_var);
    plot->Draw(); 
    plot->Write(); 


    roo_fits_file->Close(); 
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


    // Time 
    plotter->build2DHistogram("cluster time v cluster size", 160, 0, 80, 10, 0, 10); 
    
    for (int cluster_size = 1; cluster_size < 10; ++cluster_size) { 
        plotter->build1DHistogram("cluster time - cluster size " + std::to_string(cluster_size), 160, 0, 80); 
    }
    
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

    plotter->build2DHistogram("cluster energy v cluster seed energy - cuts: time, size, no edge", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size, no edge")->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->get2DHistogram("cluster energy v cluster seed energy - cuts: time, size, no edge")->GetYaxis()->SetTitle("Cluster seed energy (GeV)");

    plotter->build1DHistogram("cluster energy - no edge", 50, 0, 1.5);
    plotter->build1DHistogram("cluster time - no edge", 160, 0, 80);

    // Plots of clusters matched to tracks //
    /////////////////////////////////////////

    plotter->build2DHistogram("cluster count - matched", 47, -23, 24, 12, -6, 6);
    plotter->get2DHistogram("cluster count - matched")->GetXaxis()->SetTitle("Crystal Index - x");
    plotter->get2DHistogram("cluster count - matched")->GetYaxis()->SetTitle("Crystal Index - y");

    plotter->build1DHistogram("cluster energy - matched", 50, 0, 1.5)->GetXaxis()->SetTitle("Cluster energy (GeV)");
    plotter->build1DHistogram("cluster time - matched", 160, 0, 80)->GetXaxis()->SetTitle("Cluster time (ns)");

    // Split into top and bottom
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
    
    for (int crystal_index = 1; crystal_index < 6; ++crystal_index) { 
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

    plotter->build2DHistogram("cluster energy v p - matched - cuts: fee", 50, 0, 1.5, 50, 0, 2.0); 
    plotter->get2DHistogram("cluster energy v p - matched - cuts: fee")->GetXaxis()->SetTitle("p (GeV)"); 
    plotter->get2DHistogram("cluster energy v p - matched - cuts: fee")->GetYaxis()->SetTitle("Cluster Energy (GeV)"); 

    // Plots of tracks //
    /////////////////////

    // Plots of tracks matched to clusters
    
    // Top
    plotter->build1DHistogram("p - matched - top - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - matched - top - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - matched - top - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - matched - top - cuts: fee", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - matched - top - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    plotter->build1DHistogram("doca - matched - top", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - matched - top", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - matched - top", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - matched - top", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("tan_lambda - matched - top", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("cos(theta) - matched - top", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    // Bottom
    plotter->build1DHistogram("p - matched - bottom - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    plotter->build1DHistogram("pt - matched - bottom - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{t} [GeV]");
    plotter->build1DHistogram("px - matched - bottom - cuts: fee", 50, -0.1, 0.2)->GetXaxis()->SetTitle("p_{x} [GeV]"); 
    plotter->build1DHistogram("py - matched - bottom - cuts: fee", 50, -0.15, 0.15)->GetXaxis()->SetTitle("p_{y} [GeV]"); 
    plotter->build1DHistogram("pz - matched - bottom - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p_{z} [GeV]");

    plotter->build1DHistogram("doca - matched - bottom", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - matched - bottom", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - matched - bottom", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - matched - bottom", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("tan_lambda - matched - bottom", 100, -0.1, 0.1)->GetXaxis()->SetTitle("tan #lambda");
    plotter->build1DHistogram("cos(theta) - matched - bottom", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    // Plots of GBL tracks matched to clusters
    
    // Top
    plotter->build1DHistogram("p - matched - gbl - top - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");
    
    plotter->build1DHistogram("doca - matched - gbl - top", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - matched - gbl - top", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - matched - gbl - top", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - matched - gbl - top", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - matched - gbl - top", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    // Bottom
    plotter->build1DHistogram("p - matched - gbl - bottom - cuts: fee", 50, 0, 2.0)->GetXaxis()->SetTitle("p [GeV]");

    plotter->build1DHistogram("doca - matched - gbl - bottom", 80, -10, 10)->GetXaxis()->SetTitle("D0 [mm]");
    plotter->build1DHistogram("z0 - matched - gbl - bottom", 80, -2, 2)->GetXaxis()->SetTitle("Z0 [mm]");
    plotter->build1DHistogram("sin(phi0) - matched - gbl - bottom", 40, -0.2, 0.2)->GetXaxis()->SetTitle("sin(#phi_{0})");
    plotter->build1DHistogram("curvature - matched - gbl - bottom", 50, -0.001, 0.001)->GetXaxis()->SetTitle("#Omega");
    plotter->build1DHistogram("cos(theta) - matched - gbl - bottom", 40, -0.1, 0.1)->GetXaxis()->SetTitle("cos #theta");

    // Plots of clusters not matched to tracks
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
    if (cluster->getClusterTime() < 42 || cluster->getClusterTime() > 49.5) return false;

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
