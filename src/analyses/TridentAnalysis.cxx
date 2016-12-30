/**
 *
 * @file TridentAnalysis.cxx
 * @brief Analysis used to look at Tridents.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 29, 2015
 *
 */

#include <TridentAnalysis.h>

TridentAnalysis::TridentAnalysis()
    : class_name("TridentAnalysis") { 
    }

TridentAnalysis::~TridentAnalysis() { 
    delete ecal_utils; 
    delete matcher;
}

void TridentAnalysis::initialize() { 

    //---------------------//
    //   Event variables   //
    //---------------------//
    tuple->addVariable("event");
    tuple->addVariable("n_positrons");
    tuple->addVariable("n_tracks");
    tuple->addVariable("n_v0");
    tuple->addVariable("p_diff"); 

    //-----------------//
    //   V0 particle   //
    //-----------------//
    tuple->addVariable("invariant_mass");  
    tuple->addVariable("v0_p");
    tuple->addVariable("v_chi2");
    tuple->addVariable("vx");
    tuple->addVariable("vy");
    tuple->addVariable("vz");

    //--------------//
    //   Electron   //
    //--------------//
    tuple->addVariable("electron_chi2"); 
    tuple->addVariable("electron_ep");
    tuple->addVariable("electron_hit_n");
    tuple->addVariable("electron_has_l1");
    tuple->addVariable("electron_has_l2");
    tuple->addVariable("electron_has_l3");
    tuple->addVariable("electron_d0");
    tuple->addVariable("electron_phi0");
    tuple->addVariable("electron_omega");
    tuple->addVariable("electron_tan_lambda");
    tuple->addVariable("electron_z0");
    tuple->addVariable("electron_p");
    tuple->addVariable("electron_px"); 
    tuple->addVariable("electron_py"); 
    tuple->addVariable("electron_pz");
    tuple->addVariable("electron_time");

    tuple->addVariable("electron_cluster_energy");
    tuple->addVariable("electron_cluster_time");
    tuple->addVariable("electron_cluster_x");
    tuple->addVariable("electron_cluster_y");
    tuple->addVariable("electron_cluster_z");

    //--------------//
    //   Positron   //
    //--------------//
    tuple->addVariable("positron_chi2"); 
    tuple->addVariable("positron_ep");
    tuple->addVariable("positron_hit_n");
    tuple->addVariable("positron_has_l1");
    tuple->addVariable("positron_has_l2");
    tuple->addVariable("positron_has_l3");
    tuple->addVariable("positron_d0");
    tuple->addVariable("positron_phi0");
    tuple->addVariable("positron_omega");
    tuple->addVariable("positron_tan_lambda");
    tuple->addVariable("positron_z0");
    tuple->addVariable("positron_p"); 
    tuple->addVariable("positron_px"); 
    tuple->addVariable("positron_py"); 
    tuple->addVariable("positron_pz");
    tuple->addVariable("positron_time");

    tuple->addVariable("positron_cluster_energy");
    tuple->addVariable("positron_cluster_time");
    tuple->addVariable("positron_cluster_x");
    tuple->addVariable("positron_cluster_y");
    tuple->addVariable("positron_cluster_z");

    //---------//
    //   Top   //
    //---------//
    tuple->addVariable("top_chi2"); 
    tuple->addVariable("top_ep");
    tuple->addVariable("top_hit_n");
    tuple->addVariable("top_has_l1");
    tuple->addVariable("top_has_l2");
    tuple->addVariable("top_has_l3");
    tuple->addVariable("top_d0");
    tuple->addVariable("top_phi0");
    tuple->addVariable("top_omega");
    tuple->addVariable("top_tan_lambda");
    tuple->addVariable("top_z0");
    tuple->addVariable("top_p");
    tuple->addVariable("top_px"); 
    tuple->addVariable("top_py"); 
    tuple->addVariable("top_pz");
    tuple->addVariable("top_time");

    tuple->addVariable("top_cluster_energy");
    tuple->addVariable("top_cluster_time");
    tuple->addVariable("top_cluster_x");
    tuple->addVariable("top_cluster_y");
    tuple->addVariable("top_cluster_z");

    //---------//
    //   Bot   //
    //---------//
    tuple->addVariable("bot_chi2"); 
    tuple->addVariable("bot_ep");
    tuple->addVariable("bot_hit_n");
    tuple->addVariable("bot_has_l1");
    tuple->addVariable("bot_has_l2");
    tuple->addVariable("bot_has_l3");
    tuple->addVariable("bot_d0");
    tuple->addVariable("bot_phi0");
    tuple->addVariable("bot_omega");
    tuple->addVariable("bot_tan_lambda");
    tuple->addVariable("bot_z0");
    tuple->addVariable("bot_p");
    tuple->addVariable("bot_px"); 
    tuple->addVariable("bot_py"); 
    tuple->addVariable("bot_pz");
    tuple->addVariable("bot_time");

    tuple->addVariable("bot_cluster_energy");
    tuple->addVariable("bot_cluster_time");
    tuple->addVariable("bot_cluster_x");
    tuple->addVariable("bot_cluster_y");
    tuple->addVariable("bot_cluster_z");

    // Enable track-cluster matching plots
    matcher->enablePlots();
    matcher->useFieldMap();  

    this->bookHistograms(); 
}

void TridentAnalysis::processEvent(HpsEvent* event) { 

    ++event_counter;
    tuple->setVariableValue("event", event->getEventNumber());

    // First, check if the event contains any GBL tracks.  If it doesn't then
    // there wont be any v0 particles created from GBL tracks so the event can
    // be skipped.
    if (event->getNumberOfGblTracks() <= 1) return;
    printDebug("Event: " + std::to_string(event->getEventNumber()));
    printDebug("# Of GBL tracks: " + std::to_string(event->getNumberOfGblTracks())); 
    tuple->setVariableValue("n_tracks", event->getNumberOfGblTracks()); 

    // Find the total number of positron tracks in the event.  Only events which
    // have at least a single positron will be processed.
    std::map<GblTrack*, std::vector<HpsParticle*>> positron_map;
    std::map<GblTrack*, int> shared_hits = buildSharedHitMap(event);
    for (int track_n = 0; track_n < event->getNumberOfGblTracks(); ++track_n) { 

        printDebug("Track " + std::to_string(track_n) + " Charge: " + std::to_string(event->getGblTrack(track_n)->getCharge()));  
        // If the GBL track has a negative charge, move on to the next one.
        printDebug("Shared hits: " + std::to_string(shared_hits[event->getGblTrack(track_n)]));  
        if (event->getGblTrack(track_n)->getCharge() == -1) continue;

        // In order to keep track of multiple v0 particles created from the same
        // positron track, a mapping between a positron track and corresponding
        // v0 particles will be used.
        positron_map[event->getGblTrack(track_n)] = {};
    }
    int n_positrons{positron_map.size()}; 
    if (n_positrons == 0) return;
    ++event_has_positron;

    if (n_positrons == 1) ++event_has_single_positron;
    tuple->setVariableValue("n_positrons", n_positrons);
    printDebug("Total positrons: " + std::to_string(n_positrons)); 



    // Get the number of target constrained V0 candidates in the event.
    int n_v0{0};
    int n_particles{event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE)};

    // Loop over the collection of target contrained V0 particles and save those
    // that pass requirements placed on the Ecal clusters.  Those particles 
    // which pass the selection are added to the positron map.
    //std::vector<HpsParticle*> v0_cand; 
    for (int particle_n = 0; particle_n < n_particles; ++particle_n) { 

        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        // Only consider particles that were created from GBL tracks.
        if (particle->getType() < 32) continue;
        ++n_v0;

        if (!ecal_utils->hasGoodClusterPair(particle)) { 
            printDebug("Failed cluster selection."); 
            continue;
        }
        GblTrack* positron = static_cast<GblTrack*>(particle->getTracks()->At(1));
        positron_map[positron].push_back(particle); 
    }
    tuple->setVariableValue("n_v0", n_v0); 
    printDebug("Total v0 particles: " + std::to_string(n_v0)); 

    // Loop over the positron map and check that there are v0 particles to 
    // process.
    double total_v0{0}; 
    bool has_candidate = false;
    for (auto& candidate : positron_map) { 
        printDebug("Positron has " + std::to_string(candidate.second.size()) + " v0 candidates."); 
        if (candidate.second.size() != 0) {
            has_candidate = true;
            total_v0 += candidate.second.size(); 
        }
    }
    if (!has_candidate) return; 

    event_has_good_cluster_pair++;
    printDebug("Total v0 that have a good cluster pair: " + std::to_string(total_v0)); 
    total_v0_good_cluster_pair += total_v0;
    //if (v0_cand.size() == 0) return;

    //HpsParticle* v0{nullptr};
    std::vector<HpsParticle*> candidates; 
    double min_v0_chi2 = 1000000000000;
    int n_v0_match{0};
    for (auto& candidate : positron_map) { 
        total_v0 = 0;
        std::vector<HpsParticle*> positron_candidates; 
        for (HpsParticle* particle : candidate.second) { 
            if (!matcher->hasGoodMatch(particle)) continue;
            total_v0++; 
            n_v0_match++;
            positron_candidates.push_back(particle);  

            // Only consider the v0 particle with the smallest chi2
            //if (particle->getVertexFitChi2() > min_v0_chi2) continue;

            //min_v0_chi2 = particle->getVertexFitChi2();
            //v0 = particle;
        }

        if (positron_candidates.size() == 2) { 
            printDebug("Positron associated with " + std::to_string(positron_candidates.size()) + " particles"); 
            GblTrack* first_track  = static_cast<GblTrack*>(positron_candidates[0]->getTracks()->At(0));
            GblTrack* second_track = static_cast<GblTrack*>(positron_candidates[1]->getTracks()->At(0));
            printDebug("First Electron p: " + std::to_string(AnalysisUtils::getMagnitude(first_track->getMomentum()))); 
            printDebug("Second Electron p: " + std::to_string(AnalysisUtils::getMagnitude(second_track->getMomentum()))); 
            if (shared_hits[first_track] >= 4 
                    && shared_hits[second_track] >= 4) { 
                printDebug("Both electron tracks have a large number of shared hits."); 
                if (first_track->getChi2() < second_track->getChi2()) {
                    positron_candidates.erase(positron_candidates.begin() + 1);
                } else { 
                    positron_candidates.erase(positron_candidates.begin()); 
                }
            } else { 
                positron_candidates.clear();
            }
            printDebug("Positron now associated with " + std::to_string(positron_candidates.size()) + " particles"); 
        } else if (positron_candidates.size() > 2) { 
            printDebug("Positron has too many possibilities!");
            _dumped_positron_count++;  
            positron_candidates.clear();
        }
        if (positron_candidates.size() == 1) {
            candidates.push_back(positron_candidates[0]); 
        }
    }

    total_v0_good_track_match += n_v0_match;

    for (HpsParticle* v0 : candidates) {

        std::vector<double> p = v0->getMomentum(); 
        double v0_p = AnalysisUtils::getMagnitude(p); 
        tuple->setVariableValue("v0_p", v0_p); 
        tuple->setVariableValue("vx", v0->getVertexPosition()[0]);
        tuple->setVariableValue("vy", v0->getVertexPosition()[1]);
        tuple->setVariableValue("vz", v0->getVertexPosition()[2]);
        tuple->setVariableValue("v_chi2", v0->getVertexFitChi2()); 
        tuple->setVariableValue("invariant_mass", v0->getMass()); 

        int electron_index = 0;
        int positron_index = 1;
        SvtTrack* electron{(SvtTrack*) v0->getTracks()->At(electron_index)}; 
        SvtTrack* positron{(SvtTrack*) v0->getTracks()->At(positron_index)};

        if (positron->getCharge() == -1) { 
            electron_index = 1;
            positron_index = 0;
            electron = (SvtTrack*) v0->getTracks()->At(electron_index);
            positron = (SvtTrack*) v0->getTracks()->At(positron_index); 
        }

        int top_index = 0;
        int bot_index = 1;
        SvtTrack* top{(SvtTrack*) v0->getTracks()->At(top_index)}; 
        SvtTrack* bot{(SvtTrack*) v0->getTracks()->At(bot_index)};

        if (bot->isTopTrack()) { 
            top_index = 1;
            bot_index = 0;
            top = (SvtTrack*) v0->getTracks()->At(top_index);
            bot = (SvtTrack*) v0->getTracks()->At(bot_index); 
        }

        double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
        double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

        tuple->setVariableValue("electron_chi2", electron->getChi2());
        tuple->setVariableValue("electron_hit_n", electron->getSvtHits()->GetEntriesFast());
        tuple->setVariableValue("electron_d0", electron->getD0());
        tuple->setVariableValue("electron_phi0", electron->getPhi0());
        tuple->setVariableValue("electron_omega", electron->getOmega());
        tuple->setVariableValue("electron_tan_lambda", electron->getTanLambda());
        tuple->setVariableValue("electron_z0", electron->getZ0());
        tuple->setVariableValue("electron_p", electron_p);
        tuple->setVariableValue("electron_px", electron->getMomentum()[0]); 
        tuple->setVariableValue("electron_py", electron->getMomentum()[1]); 
        tuple->setVariableValue("electron_pz", electron->getMomentum()[2]);
        tuple->setVariableValue("electron_time", electron->getTrackTime()); 
        tuple->setVariableValue("positron_chi2", positron->getChi2());
        tuple->setVariableValue("positron_hit_n", positron->getSvtHits()->GetEntriesFast());
        tuple->setVariableValue("positron_d0", positron->getD0());
        tuple->setVariableValue("positron_phi0", positron->getPhi0());
        tuple->setVariableValue("positron_omega", positron->getOmega());
        tuple->setVariableValue("positron_tan_lambda", positron->getTanLambda());
        tuple->setVariableValue("positron_z0", positron->getZ0());
        tuple->setVariableValue("positron_p", positron_p);
        tuple->setVariableValue("positron_px", positron->getMomentum()[0]); 
        tuple->setVariableValue("positron_py", positron->getMomentum()[1]); 
        tuple->setVariableValue("positron_time", positron->getTrackTime());

        tuple->setVariableValue("p_diff", electron_p-positron_p);

        // Loop over all hits associated composing a track and check if it has a 
        // layer 1 hit.
        tuple->setVariableValue("electron_has_l1", 0);
        tuple->setVariableValue("electron_has_l2", 0);
        tuple->setVariableValue("electron_has_l3", 0);
        TRefArray* hits = electron->getSvtHits(); 
        for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
            SvtHit* hit = (SvtHit*) hits->At(hit_index); 
            if (hit->getLayer() == 1) tuple->setVariableValue("electron_has_l1", 1);
            if (hit->getLayer() == 2) tuple->setVariableValue("electron_has_l2", 1);
            if (hit->getLayer() == 3) tuple->setVariableValue("electron_has_l3", 1);
        }

        tuple->setVariableValue("positron_has_l1", 0);
        tuple->setVariableValue("positron_has_l2", 0);
        tuple->setVariableValue("positron_has_l3", 0);
        hits = positron->getSvtHits(); 
        for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
            SvtHit* hit = (SvtHit*) hits->At(hit_index); 
            if (hit->getLayer() == 1) tuple->setVariableValue("positron_has_l1", 1);
            if (hit->getLayer() == 2) tuple->setVariableValue("positron_has_l2", 1);
            if (hit->getLayer() == 3) tuple->setVariableValue("positron_has_l3", 1);
        }

        double top_p = AnalysisUtils::getMagnitude(top->getMomentum());
        double bot_p = AnalysisUtils::getMagnitude(bot->getMomentum());

        tuple->setVariableValue("top_chi2", top->getChi2());
        tuple->setVariableValue("top_hit_n", top->getSvtHits()->GetEntriesFast());
        tuple->setVariableValue("top_p", top_p);
        tuple->setVariableValue("top_d0", top->getD0());
        tuple->setVariableValue("top_phi0", top->getPhi0());
        tuple->setVariableValue("top_omega", top->getOmega());
        tuple->setVariableValue("top_tan_lambda", top->getTanLambda());
        tuple->setVariableValue("top_z0", top->getZ0());
        tuple->setVariableValue("top_px", top->getMomentum()[0]); 
        tuple->setVariableValue("top_py", top->getMomentum()[1]); 
        tuple->setVariableValue("top_pz", top->getMomentum()[2]);
        tuple->setVariableValue("top_time", top->getTrackTime()); 
        tuple->setVariableValue("bot_chi2", bot->getChi2());
        tuple->setVariableValue("bot_hit_n", bot->getSvtHits()->GetEntriesFast());
        tuple->setVariableValue("bot_d0", bot->getD0());
        tuple->setVariableValue("bot_phi0", bot->getPhi0());
        tuple->setVariableValue("bot_omega", bot->getOmega());
        tuple->setVariableValue("bot_tan_lambda", bot->getTanLambda());
        tuple->setVariableValue("bot_z0", bot->getZ0());
        tuple->setVariableValue("bot_p", bot_p);
        tuple->setVariableValue("bot_px", bot->getMomentum()[0]); 
        tuple->setVariableValue("bot_py", bot->getMomentum()[1]); 
        tuple->setVariableValue("bot_time", bot->getTrackTime());

        // Loop over all hits associated composing a track and check if it has a 
        // layer 1 hit.
        tuple->setVariableValue("top_has_l1", 0);
        tuple->setVariableValue("top_has_l2", 0);
        tuple->setVariableValue("top_has_l3", 0);
        hits = top->getSvtHits(); 
        for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
            SvtHit* hit = (SvtHit*) hits->At(hit_index); 
            if (hit->getLayer() == 1) tuple->setVariableValue("top_has_l1", 1);
            if (hit->getLayer() == 2) tuple->setVariableValue("top_has_l2", 1);
            if (hit->getLayer() == 3) tuple->setVariableValue("top_has_l3", 1);
        }

        tuple->setVariableValue("bot_has_l1", 0);
        tuple->setVariableValue("bot_has_l2", 0);
        tuple->setVariableValue("bot_has_l3", 0);
        hits = bot->getSvtHits(); 
        for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
            SvtHit* hit = (SvtHit*) hits->At(hit_index); 
            if (hit->getLayer() == 1) tuple->setVariableValue("bot_has_l1", 1);
            if (hit->getLayer() == 2) tuple->setVariableValue("bot_has_l2", 1);
            if (hit->getLayer() == 3) tuple->setVariableValue("bot_has_l3", 1);
        }


        if (v0->getClusters()->GetSize() == 2) {

            EcalCluster* electron_cluster = (EcalCluster*) v0->getClusters()->At(electron_index);
            EcalCluster* positron_cluster = (EcalCluster*) v0->getClusters()->At(positron_index);

            tuple->setVariableValue("electron_cluster_energy", electron_cluster->getEnergy());
            tuple->setVariableValue("electron_cluster_time",   electron_cluster->getClusterTime());
            tuple->setVariableValue("electron_cluster_x", electron_cluster->getPosition()[0]);
            tuple->setVariableValue("electron_cluster_y", electron_cluster->getPosition()[1]);
            tuple->setVariableValue("electron_cluster_z", electron_cluster->getPosition()[2]);

            tuple->setVariableValue("positron_cluster_energy", positron_cluster->getEnergy());
            tuple->setVariableValue("positron_cluster_time",   positron_cluster->getClusterTime());
            tuple->setVariableValue("positron_cluster_x", positron_cluster->getPosition()[0]);
            tuple->setVariableValue("positron_cluster_y", positron_cluster->getPosition()[1]);
            tuple->setVariableValue("positron_cluster_z", positron_cluster->getPosition()[2]);

            EcalCluster* top_cluster = (EcalCluster*) v0->getClusters()->At(top_index);
            EcalCluster* bot_cluster = (EcalCluster*) v0->getClusters()->At(bot_index);

            tuple->setVariableValue("top_cluster_energy", top_cluster->getEnergy());
            tuple->setVariableValue("top_cluster_time",   top_cluster->getClusterTime());
            tuple->setVariableValue("top_cluster_x", top_cluster->getPosition()[0]);
            tuple->setVariableValue("top_cluster_y", top_cluster->getPosition()[1]);
            tuple->setVariableValue("top_cluster_z", top_cluster->getPosition()[2]);

            tuple->setVariableValue("bot_cluster_energy", bot_cluster->getEnergy());
            tuple->setVariableValue("bot_cluster_time",   bot_cluster->getClusterTime());
            tuple->setVariableValue("bot_cluster_x", bot_cluster->getPosition()[0]);
            tuple->setVariableValue("bot_cluster_y", bot_cluster->getPosition()[1]);
            tuple->setVariableValue("bot_cluster_z", bot_cluster->getPosition()[2]);

        }
        tuple->fill();
    } 
}

void TridentAnalysis::finalize() {

    tuple->close(); 
    std::cout << std::fixed;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "% Total events processed: " << event_counter << std::endl;
    std::cout << "% Total events with a positron track: " << event_has_positron << std::endl;
    std::cout << "% Events with a single positron track: " << event_has_single_positron << std::endl;
    std::cout << "% Events with a good cluster pair: " << event_has_good_cluster_pair << std::endl;
    std::cout << "% Total v0 particles with a good cluster pair: " << total_v0_good_cluster_pair << std::endl;
    std::cout << "% Total v0 particles with a good track match: " << total_v0_good_track_match << std::endl;
    std::cout << "% Total positrons dumped: " << _dumped_positron_count << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
}

void TridentAnalysis::bookHistograms() {  
}

std::map<GblTrack*, int> TridentAnalysis::buildSharedHitMap(HpsEvent* event) { 
    std::map<GblTrack*, int> shared_hit_map;
    for (int first_trk_index = 0; first_trk_index < event->getNumberOfGblTracks(); ++first_trk_index) {
        printDebug("f Track " + std::to_string(first_trk_index));
        GblTrack* first_track = event->getGblTrack(first_trk_index);
        TRefArray* first_trk_hits = first_track->getSvtHits();
        shared_hit_map[first_track] = 0;
        for (int second_trk_index = 0; 
                second_trk_index < event->getNumberOfGblTracks(); ++second_trk_index) { 
            printDebug("s Track" + std::to_string(second_trk_index));
            GblTrack* second_track = event->getGblTrack(second_trk_index);
            if (first_track == second_track) continue;
            if (first_track->getTanLambda()*second_track->getTanLambda() < 0) continue;
            TRefArray* second_trk_hits = second_track->getSvtHits();
            int shared_hits = 0; 
            for (int first_hit_index = 0; first_hit_index < first_trk_hits->GetSize(); ++first_hit_index) { 
                printDebug("f Hit " + std::to_string(first_hit_index));
                SvtHit* first_hit = (SvtHit*) first_trk_hits->At(first_hit_index);
                printDebug("f Time: " + std::to_string(first_hit->getTime()));
                for (int second_hit_index = 0; second_hit_index < second_trk_hits->GetSize(); ++second_hit_index) { 
                    printDebug("s Hit " + std::to_string(second_hit_index));
                    SvtHit* second_hit = (SvtHit*) second_trk_hits->At(second_hit_index);
                    if (first_hit->getPosition()[2] == second_hit->getPosition()[2]) { 
                        printDebug("s Time: " + std::to_string(second_hit->getTime()));
                        if (first_hit->getTime() == second_hit->getTime()) shared_hits++;
                    }
                }
            }
            if(shared_hits > shared_hit_map[first_track]) shared_hit_map[first_track] = shared_hits;
        }
    }
    return shared_hit_map;
}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}

void TridentAnalysis::printDebug(std::string message) {
    if (_debug) std::cout << "[ TridentAnalysis ]: " << message << std::endl;
}
