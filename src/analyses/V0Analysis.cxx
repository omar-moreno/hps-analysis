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
      tuple(new FlatTupleMaker("trident_analysis.root", "results")), 
      matcher(new TrackClusterMatcher()),  
      class_name("V0Analysis") {
}

V0Analysis::~V0Analysis() { 
    delete ecal_utils;
    //delete tuple;
    delete matcher; 
}

void V0Analysis::initialize() {
    tuple->addVariable("cluster_energy_high");
    tuple->addVariable("cluster_energy_low");
    tuple->addVariable("cluster_x_high");
    tuple->addVariable("cluster_y_high");
    tuple->addVariable("cluster_z_high");
    tuple->addVariable("cluster_x_low");
    tuple->addVariable("cluster_y_low");
    tuple->addVariable("cluster_z_low");
    tuple->addVariable("p_sum");
    tuple->addVariable("electron_px"); 
    tuple->addVariable("electron_py"); 
    tuple->addVariable("electron_pz");
    tuple->addVariable("electron_chi2");  
    tuple->addVariable("positron_px"); 
    tuple->addVariable("positron_py"); 
    tuple->addVariable("positron_pz"); 
    tuple->addVariable("positron_chi2"); 
    tuple->addVariable("invariant_mass");  
    tuple->addVariable("vx");
    tuple->addVariable("vy");
    tuple->addVariable("vz");
    tuple->addVariable("v_chi2");
}

void V0Analysis::processEvent(HpsEvent* event) { 

    /*
     * std::cout << "[ V0Analysis ]: Event: " << event->getEventNumber() << std::endl;
     */

    // Loop over the collection of target contrained V0 particles.
    for (int particle_n = 0; particle_n < event->getNumberOfParticles(HpsParticle::TC_V0_CANDIDATE); ++particle_n) {

        // Get the nth V0 particle from the event.
        HpsParticle* particle = event->getParticle(HpsParticle::TC_V0_CANDIDATE, particle_n);

        /*
         * std::cout << "[ V0Analysis ]: Number of tracks: " << particle->getTracks()->GetEntriesFast() << std::endl;
         * std::cout << "[ V0Analysis ]: Number of clusters: " << particle->getClusters()->GetEntriesFast() << std::endl;
         */

        // Require the particle to have been created using tracks that were 
        // refit with GBL.
        if (particle->getType() < 32) {
            /*
             * std::cout << "[ V0Analysis ]: Skipping seed track particle." << std::endl;
             */
            continue; 
        }
        /*
         * std::cout << "[ V0Analysis ]: Particle type: " << particle->getType() << std::endl;
         */

        // Require each of the tracks associated with a V0 particle to have 
        // Ecal clusters matched to them.  Also require that the two clusters:
        // 1) Are coincident to within 1.6 ns
        // 2) Are in opposite Ecal volumes
        if (!ecal_utils->hasGoodClusterPair(particle)) continue;
        /*
         * std::cout << "[ V0Analysis ]: V0 particle with good Ecal cluster pair found." << std::endl; 
         */

        // Get the daughter particles
        TRefArray* daughter_particles = particle->getParticles(); 

        // Check that each of the particles has a good track match.
        if (!matcher->isGoodMatch((HpsParticle*) daughter_particles->At(0)) 
                ||!matcher->isGoodMatch((HpsParticle*) daughter_particles->At(1))) continue; 

        SvtTrack* electron = (SvtTrack*) particle->getTracks()->At(0); 
        SvtTrack* positron = (SvtTrack*) particle->getTracks()->At(1);
        if (positron->getCharge() == -1) { 
            electron = (SvtTrack*) particle->getTracks()->At(1);
            positron = (SvtTrack*) particle->getTracks()->At(0); 
        }

        EcalCluster* cluster_high = (EcalCluster*) particle->getClusters()->At(0);
        EcalCluster* cluster_low  = (EcalCluster*) particle->getClusters()->At(1);
        if (cluster_high->getEnergy() < cluster_low->getEnergy()) { 
            cluster_high = (EcalCluster*) particle->getClusters()->At(1);
            cluster_low  = (EcalCluster*) particle->getClusters()->At(0);
        }

        tuple->setVariableValue("cluster_energy_high", cluster_high->getEnergy());
        tuple->setVariableValue("cluster_energy_low", cluster_low->getEnergy());
        tuple->setVariableValue("cluster_x_high", cluster_high->getPosition()[0]);
        tuple->setVariableValue("cluster_y_high", cluster_high->getPosition()[1]);
        tuple->setVariableValue("cluster_z_high", cluster_high->getPosition()[2] );
        tuple->setVariableValue("cluster_x_low",  cluster_low->getPosition()[0]);
        tuple->setVariableValue("cluster_y_low",  cluster_low->getPosition()[1]);
        tuple->setVariableValue("cluster_z_low",  cluster_low->getPosition()[2]);

        // Calculate the momentum of the electron and positrons 
        std::vector<double> p = particle->getMomentum(); 
        double p_sum = AnalysisUtils::getMagnitude(p); 
        double electron_p = AnalysisUtils::getMagnitude(electron->getMomentum());
        double positron_p = AnalysisUtils::getMagnitude(positron->getMomentum());

        tuple->setVariableValue("p_sum", p_sum); 
        tuple->setVariableValue("electron_px", electron->getMomentum()[0]); 
        tuple->setVariableValue("electron_py", electron->getMomentum()[1]); 
        tuple->setVariableValue("electron_pz", electron->getMomentum()[2]); 
        tuple->setVariableValue("positron_px", positron->getMomentum()[0]); 
        tuple->setVariableValue("positron_py", positron->getMomentum()[1]); 
        tuple->setVariableValue("positron_pz", positron->getMomentum()[2]);
        tuple->setVariableValue("electron_chi2", electron->getChi2());
        tuple->setVariableValue("positron_chi2", positron->getChi2());
        tuple->setVariableValue("vx", particle->getVertexPosition()[0]);
        tuple->setVariableValue("vy", particle->getVertexPosition()[1]);
        tuple->setVariableValue("vz", particle->getVertexPosition()[2]);
        tuple->setVariableValue("invariant_mass", particle->getMass()); 

        tuple->fill();
    }

}

void V0Analysis::finalize() { 
    tuple->close(); 
}

void V0Analysis::bookHistograms() { 

}

std::string V0Analysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}

