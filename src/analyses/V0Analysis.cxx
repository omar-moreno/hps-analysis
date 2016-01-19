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
      class_name("V0Analysis") {
}

V0Analysis::~V0Analysis() { 
    delete ecal_utils;
}

void V0Analysis::initialize() { 
}

void V0Analysis::processEvent(HpsEvent* event) { 

    // Loop over the collection of target contrained V0 particles.
    for (int particle_n = 0; particle_n < event->getNumberOfParticles(HpsParticle::UC_V0_CANDIDATE); ++particle_n) {
        
        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::UC_V0_CANDIDATE, particle_n);

        // Require each of the tracks associated with a V0 particle to have 
        // Ecal clusters matched to them.  Also require that the two clusters:
        // 1) Are coincident to within 2 ns
        // 2) Are in opposite Ecal volumes
        if (!ecal_utils->hasGoodClusterPair(particle)) continue;
        
        std::cout << "V0 particle with good Ecal cluster pair found." << std::endl; 
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

