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
    if (n_positrons == 0) return;

    // Get the number of target constrained V0 candidates in the event.
    int n_particles = event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE);
    tuple->setVariableValue("n_v0", n_particles); 
    
    // Loop over the collection of target contrained V0 particles.
    for (int particle_n = 0; particle_n < n_particles; ++particle_n) { 
    
        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        // Only consider particles that were created from GBL tracks.
        if (particle->getType() < 32) continue;
        
        // Require each of the tracks associated with a V0 particle to have 
        // Ecal clusters matched to them.  Also require that the two clusters:
        // 1) Are coincident to within 1.6 ns
        // 2) Are in opposite Ecal volumes
        if (!ecal_utils->hasGoodClusterPair(particle)) continue;

    }
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
