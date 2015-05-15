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
      canvas(NULL), 
      output_file(new TFile("pairs_analysis_results.root", "RECREATE")),
      class_name("PairsAnalysis") { 
}

PairsAnalysis::~PairsAnalysis() { 

    delete canvas; 
    delete output_file;

    uc_vtx_plots.clear();
}

void PairsAnalysis::initialize() { 

    this->bookHistograms(); 
}

void PairsAnalysis::processEvent(HpsEvent* event) { 

    //if (event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES) == 0) return;

    uc_vtx_plots["Number of UC VTX particles"]->Fill(event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES));

    for (int particle_n = 0; 
            particle_n < event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES);  
            ++particle_n) { 

        particle = event->getParticle(HpsEvent::UC_VTX_PARTICLES, particle_n);

        uc_vtx_plots["UC VTX invariant mass"]->Fill(particle->Mass());
        uc_vtx_plots["UC VTX Vx"]->Fill(particle->getVertexPosition()[0]); 
        uc_vtx_plots["UC VTX Vy"]->Fill(particle->getVertexPosition()[1]); 
        uc_vtx_plots["UC VTX Vz"]->Fill(particle->getVertexPosition()[2]); 

    }
}

void PairsAnalysis::finalize() { 

    canvas->Print("pairs_analysis_results.pdf[");
    std::unordered_map<std::string, TH1F*>::iterator uc_vtx_plots_it = uc_vtx_plots.begin();
    for (uc_vtx_plots_it; uc_vtx_plots_it != uc_vtx_plots.end(); ++uc_vtx_plots_it) { 
        uc_vtx_plots_it->second->Draw();
        uc_vtx_plots_it->second->Write();
        canvas->Print("pairs_analysis_results.pdf(");
    }
    canvas->Print("pairs_analysis_results.pdf]");
    output_file->Close();
}

void PairsAnalysis::bookHistograms() { 

    canvas = new TCanvas("canvas", "canvas", 500, 500);

    uc_vtx_plots["Number of UC VTX particles"] 
        = new TH1F("n_uc_vtx_particles", "n_uc_vtx_particles", 5, 0, 5);
    uc_vtx_plots["UC VTX invariant mass"] 
        = new TH1F("uc_vtx_invariant_mass", "uc_vtx_invariant_mass", 100, 0, 0.200);
    uc_vtx_plots["UC VTX Vx"] = new TH1F("uc_vtx_vx", "uc_vtx_vx", 50, -1, 1); 
    uc_vtx_plots["UC VTX Vy"] = new TH1F("uc_vtx_vy", "uc_vtx_vy", 50, -0.6, 0.6); 
    uc_vtx_plots["UC VTX Vz"] = new TH1F("uc_vtx_vz", "uc_vtx_vz", 50, -10, 10);
    uc_vtx_plots["UC VTX chi2"] = new TH1F("uc_vtx_chi2", "uc_vtx_chi2", 25, 0, 25); 
}

std::string PairsAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
