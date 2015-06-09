/**
 * @file PairsAnalysis.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 14, 2015
 *
 */

#include <PairsAnalysis.h>

PairsAnalysis::PairsAnalysis()
    : particle(NULL),
      uc_vtx_plotter(new Plotter()),
      uc_vtx_epem_plotter(new Plotter()),
      uc_vtx_emem_plotter(new Plotter()),
      class_name("PairsAnalysis") { 
}

PairsAnalysis::~PairsAnalysis() { 
    delete uc_vtx_epem_plotter;
    delete uc_vtx_emem_plotter;
}

void PairsAnalysis::initialize() { 
    this->bookHistograms(); 
}

void PairsAnalysis::processEvent(HpsEvent* event) { 

    //if (!event->isPair1Trigger()) return;

    if (event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES) == 0) return;

    uc_vtx_plotter->get1DHistogram("number of vtx particles")->Fill
        (event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES));

    for (int particle_n = 0; 
            particle_n < event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES);  
            ++particle_n) { 

        particle = event->getParticle(HpsEvent::UC_VTX_PARTICLES, particle_n);


        TRefArray* daughter_particles = particle->getParticles();
        int charge_sum = 0;
        int charge[2] = {0, 0};
        double p_mag[2] = {0, 0};
        for (int d_particle_n = 0; d_particle_n <daughter_particles->GetSize(); ++d_particle_n) { 
            charge_sum += ((HpsParticle*) daughter_particles->At(d_particle_n))->getCharge();
           
            std::vector<double> p
                = ((HpsParticle*) daughter_particles->At(d_particle_n))->getMomentum();

            if (daughter_particles->GetSize() == 2) { 
               
                    charge[d_particle_n] 
                        = ((HpsParticle*) daughter_particles->At(d_particle_n))->getCharge();
                    p_mag[d_particle_n] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
            }
        }
    
        std::vector<double> p = particle->getMomentum();
        double v0_p_mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double v0_pt = sqrt(p[0]*p[0] + p[1]*p[1]);
        
        if (charge_sum == 0 && daughter_particles->GetSize() == 2 && v0_p_mag < 1.1*1.056 &&
                p[0] < .2 && p[1] < .2 && abs(particle->getVertexPosition()[0]) < 2 && abs(particle->getVertexPosition()[1]) < 1.0 
                && abs(particle->getVertexPosition()[3]) < 25) { 
            
            uc_vtx_epem_plotter->get1DHistogram("invariant mass")->Fill(particle->Mass());
            uc_vtx_epem_plotter->get1DHistogram("vx")->Fill(particle->getVertexPosition()[0]); 
            uc_vtx_epem_plotter->get1DHistogram("vy")->Fill(particle->getVertexPosition()[1]); 
            uc_vtx_epem_plotter->get1DHistogram("vz")->Fill(particle->getVertexPosition()[2]); 

            if (charge[0] < 0 && charge[1] > 0) { 
                uc_vtx_epem_plotter->get2DHistogram("p[e+] v p[e-]")->Fill(p_mag[1], p_mag[0]); 
            } else if (charge[0] > 0 && charge[1] < 0) { 
                uc_vtx_epem_plotter->get2DHistogram("p[e+] v p[e-]")->Fill(p_mag[0], p_mag[1]); 
            }

            uc_vtx_epem_plotter->get1DHistogram("p")->Fill(v0_p_mag);
            uc_vtx_epem_plotter->get1DHistogram("pt")->Fill(v0_pt);
            uc_vtx_epem_plotter->get1DHistogram("px")->Fill(p[0]);
            uc_vtx_epem_plotter->get1DHistogram("py")->Fill(p[1]);
            uc_vtx_epem_plotter->get1DHistogram("pz")->Fill(p[2]);

        } else if (daughter_particles->GetSize() == 2 && charge[0] < 0 && charge[1] < 0) { 
            
            uc_vtx_emem_plotter->get1DHistogram("invariant mass")->Fill(particle->Mass());
            uc_vtx_emem_plotter->get1DHistogram("vx")->Fill(particle->getVertexPosition()[0]); 
            uc_vtx_emem_plotter->get1DHistogram("vy")->Fill(particle->getVertexPosition()[1]); 
            uc_vtx_emem_plotter->get1DHistogram("vz")->Fill(particle->getVertexPosition()[2]); 
                
            uc_vtx_emem_plotter->get2DHistogram("p[e-] v p[e-]")->Fill(p_mag[1], p_mag[0]);
        } 
    }
}

void PairsAnalysis::finalize() { 

    uc_vtx_plotter->saveToRootFile("all_pairs_analysis.root");
    uc_vtx_plotter->saveToPdf("all_pairs_analysis.pdf");
    
    uc_vtx_epem_plotter->saveToRootFile("epem_analysis.root");
    uc_vtx_epem_plotter->saveToPdf("epem_analysis.pdf");
    
    uc_vtx_emem_plotter->saveToRootFile("emem_analysis.root");
    uc_vtx_emem_plotter->saveToPdf("emem_analysis.pdf");
}

void PairsAnalysis::bookHistograms() { 

    uc_vtx_plotter->setType("float");
    uc_vtx_plotter->build1DHistogram("number of vtx particles", 10, 0, 10);

    uc_vtx_epem_plotter->setType("float");
    uc_vtx_epem_plotter->build1DHistogram("invariant mass", 100, 0, 0.200);
    uc_vtx_epem_plotter->build1DHistogram("vx", 50, -3, 3); 
    uc_vtx_epem_plotter->build1DHistogram("vy", 50, -2, 2); 
    uc_vtx_epem_plotter->build1DHistogram("vz", 80, -40, 40);
    uc_vtx_epem_plotter->build1DHistogram("#chi^{2}", 25, 0, 25); 
    uc_vtx_epem_plotter->build1DHistogram("p", 50, 0, 2.0);
    uc_vtx_epem_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    uc_vtx_epem_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    uc_vtx_epem_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    uc_vtx_epem_plotter->build1DHistogram("pz", 50, 0, 2.0);

    uc_vtx_emem_plotter->setType("float");
    uc_vtx_emem_plotter->build1DHistogram("invariant mass", 100, 0, 0.200);
    uc_vtx_emem_plotter->build1DHistogram("vx", 50, -3, 3); 
    uc_vtx_emem_plotter->build1DHistogram("vy", 50, -2, 2); 
    uc_vtx_emem_plotter->build1DHistogram("vz", 80, -40, 40);
    uc_vtx_emem_plotter->build1DHistogram("#chi^{2}", 25, 0, 25); 
    uc_vtx_emem_plotter->build1DHistogram("p", 50, 0, 2.0);
    uc_vtx_emem_plotter->build1DHistogram("pt", 50, -0.1, 0.2);
    uc_vtx_emem_plotter->build1DHistogram("px", 50, -0.1, 0.2); 
    uc_vtx_emem_plotter->build1DHistogram("py", 50, -0.15, 0.15); 
    uc_vtx_emem_plotter->build1DHistogram("pz", 50, 0, 2.0);

    uc_vtx_epem_plotter->build2DHistogram("p[e+] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
    uc_vtx_emem_plotter->build2DHistogram("p[e-] v p[e-]", 50, 0, 2.0, 50, 0, 2.0);
}

std::string PairsAnalysis::toString() { 
    
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
