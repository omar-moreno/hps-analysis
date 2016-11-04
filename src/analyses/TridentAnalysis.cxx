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

    // Enable track-cluster matching plots
    matcher->enablePlots();
    matcher->useFieldMap();  

    this->bookHistograms(); 
}

void TridentAnalysis::processEvent(HpsEvent* event) { 

    ++event_counter;
    tuple->setVariableValue("event", event->getEventNumber());
   
    // If the event doesn't have any tracks, skip the rest of the processing 
    // chain.
    double n_tracks = event->getNumberOfGblTracks();
    if (n_tracks == 0) return;
    ++event_has_track;
    tuple->setVariableValue("n_tracks", n_tracks); 

    // Get the number of positrons in the event.  If the event doesn't have
    // any positron tracks, skip the rest of the processing chain.
    double n_positrons = 0;
    for (int track_n = 0; track_n < event->getNumberOfGblTracks(); ++track_n) { 
        if (event->getGblTrack(track_n)->getCharge() == 1) n_positrons++;
    }
    if (n_positrons == 0) return;
    ++event_has_positron;

    if (n_positrons != 1) return;
    ++event_has_single_positron;
    tuple->setVariableValue("n_positrons", n_positrons);

    // Get the number of target constrained V0 candidates in the event.
    int n_v0{0};
    int n_particles{event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE)};

    // Loop over the collection of target contrained V0 particles.
    HpsParticle* v0{nullptr};
    double min_v0_chi2 = 1000000000000;
    for (int particle_n = 0; particle_n < n_particles; ++particle_n) { 
    
        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        // Only consider particles that were created from GBL tracks.
        if (particle->getType() < 32) continue;
        ++n_v0; 

        // Require the two tracks associated with the v0 particle to be 
        // oppositely charged.      
        SvtTrack* first_track{(SvtTrack*) particle->getTracks()->At(0)};
        SvtTrack* second_track{(SvtTrack*) particle->getTracks()->At(1)};
        if (first_track->getCharge()*second_track->getCharge() == 1) continue;
        
        // Require the two tracks associated with the v0 particle to be in 
        // opposite volumes.
        if (first_track->getMomentum()[1]*second_track->getMomentum()[1] > 0) continue;

        // Only consider the v0 particle with the smallest chi2
        if (particle->getVertexFitChi2() > min_v0_chi2) continue;

        min_v0_chi2 = particle->getVertexFitChi2();
        v0 = particle; 
    }
    tuple->setVariableValue("n_v0", n_v0); 
    
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
    tuple->setVariableValue("electron_hit_n", electron->getSvtHits()->GetEntriesFast());
    tuple->setVariableValue("electron_p", electron_p);
    tuple->setVariableValue("electron_px", electron->getMomentum()[0]); 
    tuple->setVariableValue("electron_py", electron->getMomentum()[1]); 
    tuple->setVariableValue("electron_pz", electron->getMomentum()[2]);
    tuple->setVariableValue("electron_time", electron->getTrackTime()); 
    tuple->setVariableValue("positron_hit_n", positron->getSvtHits()->GetEntriesFast());
    tuple->setVariableValue("positron_p", positron_p);
    tuple->setVariableValue("positron_px", positron->getMomentum()[0]); 
    tuple->setVariableValue("positron_py", positron->getMomentum()[1]); 
    tuple->setVariableValue("positron_time", positron->getTrackTime());

    // Loop over all hits associated composing a track and check if it has a 
    // layer 1 hit.
    tuple->setVariableValue("electron_has_l1", 0);
    TRefArray* hits = electron->getSvtHits(); 
    for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
        SvtHit* hit = (SvtHit*) hits->At(hit_index); 
        if (hit->getLayer() == 1) tuple->setVariableValue("electron_has_l1", 1);
    }

    tuple->setVariableValue("positron_has_l1", 0);
    hits = positron->getSvtHits(); 
    for (int hit_index = 0; hit_index < hits->GetEntriesFast(); ++hit_index) { 
        SvtHit* hit = (SvtHit*) hits->At(hit_index); 
        if (hit->getLayer() == 1) tuple->setVariableValue("positron_has_l1", 1);
    }

    if (v0->getClusters()->GetSize() == 2) {
    
        EcalCluster* electron_cluster = (EcalCluster*) v0->getClusters()->At(electron_index);
        EcalCluster* positron_cluster = (EcalCluster*) v0->getClusters()->At(positron_index);
    
        tuple->setVariableValue("electron_cluster_energy", electron_cluster->getEnergy());
        tuple->setVariableValue("electron_cluster_time",   electron_cluster->getClusterTime());
        tuple->setVariableValue("electron_cluster_x", electron_cluster->getPosition()[0]);
        tuple->setVariableValue("electron_cluster_y", electron_cluster->getPosition()[1]);
        tuple->setVariableValue("electron_cluster_z", electron_cluster->getPosition()[2] );

        tuple->setVariableValue("positron_cluster_energy", positron_cluster->getEnergy());
        tuple->setVariableValue("positron_cluster_time",   positron_cluster->getClusterTime());
        tuple->setVariableValue("positron_cluster_x", positron_cluster->getPosition()[0]);
        tuple->setVariableValue("positron_cluster_y", positron_cluster->getPosition()[1]);
        tuple->setVariableValue("positron_cluster_z", positron_cluster->getPosition()[2] );
    } 
    tuple->fill();
}

void TridentAnalysis::finalize() {
    
    tuple->close(); 
    std::cout << std::fixed;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "% Total events processed: " << event_counter << std::endl;
    std::cout << "% Total events with a track: " << event_has_track << std::endl;
    std::cout << "% Total events with a positron track: " << event_has_positron << std::endl;
    std::cout << "% Events with a single positron track: " << event_has_single_positron << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
}

void TridentAnalysis::bookHistograms() {  
}

std::string TridentAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
