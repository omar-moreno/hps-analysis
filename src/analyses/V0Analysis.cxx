/**
 *
 * @file V0Analysis.h
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date January 18, 2016
 * @brief Analysis that uses V0 particles to look at Tridents.
 *
 */

#include <V0Analysis.h>

V0Analysis::V0Analysis() 
    : ecal_utils(new EcalUtils()),
      matcher(new TrackClusterMatcher()),  
      class_name("V0Analysis") {
}

V0Analysis::~V0Analysis() { 
    delete ecal_utils;
}

void V0Analysis::initialize() { 
}

void V0Analysis::processEvent(HpsEvent* event) { 

    std::cout << "[ V0Analysis ]: Event: " << event->getEventNumber() << std::endl;

    // Loop over the collection of target contrained V0 particles.
    for (int particle_n = 0; particle_n < event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE); ++particle_n) {

        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        std::cout << "[ V0Analysis ]: Number of tracks: " << particle->getTracks()->GetEntriesFast() << std::endl;
        std::cout << "[ V0Analysis ]: Number of clusters: " << particle->getClusters()->GetEntriesFast() << std::endl;

        // Require the particle to have been created using tracks that were 
        // refit with GBL.
        if (particle->getType() < 32) {
            std::cout << "[ V0Analysis ]: Skipping seed track particle." << std::endl;
            continue; 
        }
        std::cout << "[ V0Analysis ]: Particle type: " << particle->getType() << std::endl;

        // Require each of the tracks associated with a V0 particle to have 
        // Ecal clusters matched to them.  Also require that the two clusters:
        // 1) Are coincident to within 1.6 ns
        // 2) Are in opposite Ecal volumes
        if (!ecal_utils->hasGoodClusterPair(particle)) continue;
        std::cout << "[ V0Analysis ]: V0 particle with good Ecal cluster pair found." << std::endl; 
    
        // Get the daughter particles
        TRefArray* daughter_particles = particle->getParticles(); 

        // Check that each of the particles has a good track match.
        if (!matcher->isGoodMatch((HpsParticle*) daughter_particles->At(0)) 
                ||!matcher->isGoodMatch((HpsParticle*) daughter_particles->At(1))) continue; 

        std::cout << "[ V0Analysis ]: V0 particle with good track-cluster match found." << std::endl; 
    }
}

void V0Analysis::finalize() { 

}

void V0Analysis::bookHistograms() { 

}

std::string V0Analysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}

