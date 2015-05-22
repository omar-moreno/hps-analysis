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

    //if (!event->isPair1Trigger()) return;

    if (event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES) == 0) return;

    uc_vtx_plots["Number of UC VTX particles"]->Fill(event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES));

    for (int particle_n = 0; 
            particle_n < event->getNumberOfParticles(HpsEvent::UC_VTX_PARTICLES);  
            ++particle_n) { 

        particle = event->getParticle(HpsEvent::UC_VTX_PARTICLES, particle_n);

        uc_vtx_plots["UC VTX invariant mass"]->Fill(particle->Mass());
        uc_vtx_plots["UC VTX Vx"]->Fill(particle->getVertexPosition()[0]); 
        uc_vtx_plots["UC VTX Vy"]->Fill(particle->getVertexPosition()[1]); 
        uc_vtx_plots["UC VTX Vz"]->Fill(particle->getVertexPosition()[2]); 
        uc_vtx_plots["UC VTX chi2"]->Fill(particle->getVertexFitChi2()); 

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
            
            uc_vtx_epem_plots["UC VTX invariant mass"]->Fill(particle->Mass());
            uc_vtx_epem_plots["UC VTX Vx"]->Fill(particle->getVertexPosition()[0]); 
            uc_vtx_epem_plots["UC VTX Vy"]->Fill(particle->getVertexPosition()[1]); 
            uc_vtx_epem_plots["UC VTX Vz"]->Fill(particle->getVertexPosition()[2]); 

            if (charge[0] < 0 && charge[1] > 0) { 
                uc_vtx_epem_plots_2d["p[e+] v p[e-]"]->Fill(p_mag[1], p_mag[0]); 
            } else if (charge[0] > 0 && charge[1] < 0) { 
                uc_vtx_epem_plots_2d["p[e+] v p[e-]"]->Fill(p_mag[0], p_mag[1]); 
            }

            uc_vtx_epem_plots["p"]->Fill(v0_p_mag);
            uc_vtx_epem_plots["pt"]->Fill(v0_pt);
            uc_vtx_epem_plots["px"]->Fill(p[0]);
            uc_vtx_epem_plots["py"]->Fill(p[1]);
            uc_vtx_epem_plots["pz"]->Fill(p[2]);

        } else if (daughter_particles->GetSize() == 2 && charge[0] < 0 && charge[1] < 0) { 
            
            uc_vtx_emem_plots["UC VTX invariant mass"]->Fill(particle->Mass());
            uc_vtx_emem_plots["UC VTX Vx"]->Fill(particle->getVertexPosition()[0]); 
            uc_vtx_emem_plots["UC VTX Vy"]->Fill(particle->getVertexPosition()[1]); 
            uc_vtx_emem_plots["UC VTX Vz"]->Fill(particle->getVertexPosition()[2]); 
                
            uc_vtx_emem_plots_2d["p[e-] v p[e-]"]->Fill(p_mag[1], p_mag[0]);
        } 
    }
}

void PairsAnalysis::finalize() { 

    canvas->Print("pairs_analysis_results.pdf[");
    std::unordered_map<std::string, TH1F*>::iterator uc_vtx_plots_it = uc_vtx_plots.begin();
    for (uc_vtx_plots_it; uc_vtx_plots_it != uc_vtx_plots.end(); ++uc_vtx_plots_it) { 
        uc_vtx_plots_it->second->Draw();
        uc_vtx_plots_it->second->Write();

        if (uc_vtx_epem_plots[uc_vtx_plots_it->first] != NULL) {
            uc_vtx_epem_plots[uc_vtx_plots_it->first]->Draw("same");
            uc_vtx_epem_plots[uc_vtx_plots_it->first]->Write();
        }

        if (uc_vtx_emem_plots[uc_vtx_plots_it->first] != NULL) { 
            uc_vtx_emem_plots[uc_vtx_plots_it->first]->Draw("same");
            uc_vtx_emem_plots[uc_vtx_plots_it->first]->Write();
        }

        canvas->Print("pairs_analysis_results.pdf(");
    }

    std::unordered_map<std::string, TH2F*>::iterator uc_vtx_epem_plots_2d_it 
        = uc_vtx_epem_plots_2d.begin();
    for (uc_vtx_epem_plots_2d_it; 
            uc_vtx_epem_plots_2d_it != uc_vtx_epem_plots_2d.end();
            ++uc_vtx_epem_plots_2d_it) { 

        uc_vtx_epem_plots_2d_it->second->Draw("colz");
        uc_vtx_epem_plots_2d_it->second->Write();
        canvas->Print("pairs_analysis_results.pdf(");
    }
   

    uc_vtx_epem_plots["p"]->Draw();
        canvas->Print("pairs_analysis_results.pdf(");
    uc_vtx_epem_plots["pt"]->Draw();
        canvas->Print("pairs_analysis_results.pdf(");
    canvas->Print("pairs_analysis_results.pdf(");
    uc_vtx_epem_plots["px"]->Draw();
        canvas->Print("pairs_analysis_results.pdf(");
    uc_vtx_epem_plots["py"]->Draw();
        canvas->Print("pairs_analysis_results.pdf(");
    uc_vtx_epem_plots["pz"]->Draw();
        canvas->Print("pairs_analysis_results.pdf(");

    uc_vtx_epem_plots["p"]->Write();
    uc_vtx_epem_plots["pt"]->Write();
    uc_vtx_epem_plots["px"]->Write();
    uc_vtx_epem_plots["py"]->Write();
    uc_vtx_epem_plots["pz"]->Write();


    std::unordered_map<std::string, TH2F*>::iterator uc_vtx_emem_plots_2d_it 
        = uc_vtx_emem_plots_2d.begin();
    for (uc_vtx_emem_plots_2d_it; 
            uc_vtx_emem_plots_2d_it != uc_vtx_emem_plots_2d.end();
            ++uc_vtx_emem_plots_2d_it) { 

        uc_vtx_emem_plots_2d_it->second->Draw("colz");
        uc_vtx_emem_plots_2d_it->second->Write();
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
    uc_vtx_plots["UC VTX Vx"] = new TH1F("uc_vtx_vx", "uc_vtx_vx", 50, -3, 3); 
    uc_vtx_plots["UC VTX Vy"] = new TH1F("uc_vtx_vy", "uc_vtx_vy", 50, -2, 2); 
    uc_vtx_plots["UC VTX Vz"] = new TH1F("uc_vtx_vz", "uc_vtx_vz", 80, -40, 40);
    uc_vtx_plots["UC VTX chi2"] = new TH1F("uc_vtx_chi2", "uc_vtx_chi2", 25, 0, 25); 

    uc_vtx_epem_plots["UC VTX invariant mass"] 
        = new TH1F("uc_vtx_invariant_mass_epem", "uc_vtx_invariant_mass_epem", 100, 0, 0.200);
    uc_vtx_epem_plots["UC VTX Vx"] = new TH1F("uc_vtx_vx_epem", "uc_vtx_vx_epem", 50, -3, 3); 
    uc_vtx_epem_plots["UC VTX Vy"] = new TH1F("uc_vtx_vy_epem", "uc_vtx_vy_epem", 50, -2, 2); 
    uc_vtx_epem_plots["UC VTX Vz"] = new TH1F("uc_vtx_vz_epem", "uc_vtx_vz_epem", 80, -40, 40);
    uc_vtx_epem_plots["UC VTX chi2"] = new TH1F("uc_vtx_chi2_epem", "uc_vtx_chi2_epem", 25, 0, 25); 
    uc_vtx_epem_plots["p"] = new TH1F("p", "p", 50, 0, 2.0);
    uc_vtx_epem_plots["pt"] = new TH1F("pt", "pt", 50, -0.1, 0.2);
    uc_vtx_epem_plots["px"] = new TH1F("px", "px", 50, -0.1, 0.2); 
    uc_vtx_epem_plots["py"] = new TH1F("py", "py", 50, -0.15, 0.15); 
    uc_vtx_epem_plots["pz"] = new TH1F("pz", "pz", 50, 0, 2.0);



    uc_vtx_emem_plots["UC VTX invariant mass"] 
        = new TH1F("uc_vtx_invariant_mass_emem", "uc_vtx_invariant_mass_emem", 100, 0, 0.200);
    uc_vtx_emem_plots["UC VTX Vx"] = new TH1F("uc_vtx_vx_emem", "uc_vtx_vx_emem", 50, -3, 3); 
    uc_vtx_emem_plots["UC VTX Vy"] = new TH1F("uc_vtx_vy_emem", "uc_vtx_vy_emem", 50, -2, 2); 
    uc_vtx_emem_plots["UC VTX Vz"] = new TH1F("uc_vtx_vz_emem", "uc_vtx_vz_emem", 80, -40, 40);
    uc_vtx_emem_plots["UC VTX chi2"] = new TH1F("uc_vtx_chi2_emem", "uc_vtx_chi2_emem", 25, 0, 25); 

    uc_vtx_epem_plots_2d["p[e+] v p[e-]"] 
        = new TH2F("pep_v_pem", "pep_v_pem", 50, 0, 2.0, 50, 0, 2.0);
    
    uc_vtx_emem_plots_2d["p[e-] v p[e-]"] 
        = new TH2F("pem_v_pem", "pem_v_pem", 50, 0, 2.0, 50, 0, 2.0);

}

std::string PairsAnalysis::toString() { 

    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer; 
}
