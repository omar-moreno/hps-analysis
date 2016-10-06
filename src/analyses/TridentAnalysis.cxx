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
    : 
      ecal_utils(new EcalUtils()), 
      matcher(new TrackClusterMatcher()), 
      class_name("TridentAnalysis"), 
      good_cluster_pair_counter(0),  
      matched_event_counter(0), 
      v0_cand_counter(0) { 
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
    tuple->addVariable("electron_p");
    tuple->addVariable("electron_px"); 
    tuple->addVariable("electron_py"); 
    tuple->addVariable("electron_pz");
    tuple->addVariable("electron_time");

    tuple->addVariable("electron_cluster_energy");
    tuple->addVariable("electron_cluster_time");
    tuple->addVariable("electron_cluster_x_high");
    tuple->addVariable("electron_cluster_y_high");
    tuple->addVariable("electron_cluster_z_high");

    //--------------//
    //   Positron   //
    //--------------//
    tuple->addVariable("positron_chi2"); 
    tuple->addVariable("positron_ep");
    tuple->addVariable("positron_hit_n");
    tuple->addVariable("positron_has_l1");
    tuple->addVariable("positron_p"); 
    tuple->addVariable("positron_px"); 
    tuple->addVariable("positron_py"); 
    tuple->addVariable("positron_pz");
    tuple->addVariable("positron_time");

    tuple->addVariable("positron_cluster_energy");
    tuple->addVariable("positron_cluster_time");
    tuple->addVariable("positron_cluster_x_high");
    tuple->addVariable("positron_cluster_y_high");
    tuple->addVariable("positron_cluster_z_high");

    // Enable track-cluster matching plots
    matcher->enablePlots();
    matcher->useFieldMap();  

    this->bookHistograms(); 
}

void TridentAnalysis::processEvent(HpsEvent* event) { 

    ++event_counter;
    tuple->setVariableValue("event", event->getEventNumber());
   
    double n_tracks = event->getNumberOfGblTracks();
    tuple->setVariableValue("n_tracks", n_tracks);

    // Get the number of positrons in the event.
    double n_positrons = 0;
    for (int track_n = 0; track_n < event->getNumberOfGblTracks(); ++track_n) { 
        if (event->getGblTrack(track_n)->getCharge() == 1) n_positrons++;
    }
    tuple->setVariableValue("n_positrons", n_positrons);

    // If the event doesn't have any positrons, then no tridents have been
    // found. Skip the event.
    if (n_positrons == 0) { 
        tuple->fill(); 
        return;
    }

    // Get the number of target constrained V0 candidates in the event.
    int n_particles = event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE);
    tuple->setVariableValue("n_v0", n_particles); 
    
    // Loop over the collection of target contrained V0 particles.
    HpsParticle* v0{nullptr};
    double min_v0_chi2 = 1000000000000;
    for (int particle_n = 0; particle_n < n_particles; ++particle_n) { 
    
        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        // Only consider particles that were created from GBL tracks.
        if (particle->getType() < 32) continue;
        
        // Require the two tracks associated with the v0 particle to be 
        // oppositely charged.      
        SvtTrack* first_track{(SvtTrack*) particle->getTracks()->At(0)};
        SvtTrack* second_track{(SvtTrack*) particle->getTracks()->At(1)};
        if (first_track->getCharge()*second_track->getCharge() == 1) continue;
        
        // Require the two tracks associated with the v0 particle to be in 
        // opposite volumes.
        //if (first_track

        // Only consider the v0 particle with the smallest chi2
        if (particle->getVertexFitChi2() > min_v0_chi2) continue;

        min_v0_chi2 = particle->getVertexFitChi2();
        v0 = particle; 
    }
    
    // If a v0 particle wasn't found, don't continue processing the event.
    if (v0 == nullptr) {
        tuple->fill();
        return;
    }

    // Calculate the momentum of the electron and positrons 
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
    
    double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
    double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

    tuple->setVariableValue("electron_chi2", electron->getChi2());
    //tuple->setVariableValue("electron_ep", electron_ep);
    tuple->setVariableValue("electron_hit_n", electron->getSvtHits()->GetEntriesFast());
    tuple->setVariableValue("electron_p", electron_p);
    tuple->setVariableValue("electron_px", electron->getMomentum()[0]); 
    tuple->setVariableValue("electron_py", electron->getMomentum()[1]); 
    tuple->setVariableValue("electron_pz", electron->getMomentum()[2]);
    tuple->setVariableValue("electron_time", electron->getTrackTime()); 
    //tuple->setVariableValue("positron_ep", positron_ep);
    tuple->setVariableValue("positron_hit_n", positron->getSvtHits()->GetEntriesFast());
    tuple->setVariableValue("positron_p", positron_p);
    tuple->setVariableValue("positron_px", positron->getMomentum()[0]); 
    tuple->setVariableValue("positron_py", positron->getMomentum()[1]); 
    tuple->setVariableValue("positron_time", positron->getTrackTime());

    if (v0->getClusters()->GetSize() == 2) {
    
        EcalCluster* electron_cluster = (EcalCluster*) v0->getClusters()->At(electron_index);
        EcalCluster* positron_cluster = (EcalCluster*) v0->getClusters()->At(positron_index);
    
        tuple->setVariableValue("electron_cluster_energy", electron_cluster->getEnergy());
        tuple->setVariableValue("electron_cluster_time",   electron_cluster->getClusterTime());
        tuple->setVariableValue("electron_cluster_x_high", electron_cluster->getPosition()[0]);
        tuple->setVariableValue("electron_cluster_y_high", electron_cluster->getPosition()[1]);
        tuple->setVariableValue("electron_cluster_z_high", electron_cluster->getPosition()[2] );

        tuple->setVariableValue("positron_cluster_energy", positron_cluster->getEnergy());
        tuple->setVariableValue("positron_cluster_time",   positron_cluster->getClusterTime());
        tuple->setVariableValue("positron_cluster_x_high", positron_cluster->getPosition()[0]);
        tuple->setVariableValue("positron_cluster_y_high", positron_cluster->getPosition()[1]);
        tuple->setVariableValue("positron_cluster_z_high", positron_cluster->getPosition()[2] );
    } 
    tuple->fill();
}

void TridentAnalysis::finalize() {
    
    // Save the track-Ecal cluster matching plots to a ROOT file
    matcher->saveHistograms();  

    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "Number of pair1 triggers: " << event_counter << std::endl;
    std::cout << "Good clusters: " << good_cluster_pair_counter << "/" <<  event_counter << " = "
              << good_cluster_pair_counter/event_counter << " %" << std::endl;
    std::cout << "Trident's passing cuts: " << v0_cand_counter << "/" << event_counter << " = " 
              << v0_cand_counter/event_counter << " %" << std::endl; 
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
}

void TridentAnalysis::bookHistograms() {  
}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
