
//------------------//
//--- C++ StdLib ---//
//------------------//
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <string>
#include <list>

//------------//
//--- ROOT ---//
//------------//
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

using namespace std; 

int main(int argc, char **argv) {

    string file_list_name;
    int event_count = -1; 

    // Parse all the command line arguments.  If there are no valid command
    // line arguments passed, print the usage and exit the application
    static struct option long_options[] = {
        {"file_list",  required_argument, 0, 'l'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int option_char; 
    while ((option_char = getopt_long(argc, argv, "l:", long_options, &option_index)) != -1) {

        switch(option_char){
            case 'l':
                file_list_name = optarg;
                break; 
            default: 
                return EXIT_FAILURE;
        }
    }

    // If a list of DST files was not specified by the user, warn the user
    // and exit the application.
    if (file_list_name.empty()) { 
        cerr << "\n[ HPS ANALYZER ]: Please specify a list of DST files to process." << endl;
        return EXIT_FAILURE;
    }

    // Create a list of files to process
    list<string> files; 
    string file;

    ifstream file_list(file_list_name.c_str(), ifstream::in);
    if (!file_list.is_open()) { 
        cerr << "\n[ HPS ANALYZER ]: Failed to open file " << file_list_name << endl;
        return EXIT_FAILURE;
    }

    while (file_list >> file) { 
        files.push_back(file); 
    }
    file_list.close();

    std::vector<TFile*> t_files; 
    std::vector<TGraph*> chi2_graphs; 
    std::vector<TGraph*> omega_graphs; 
    std::vector<TGraph*> phi_graphs; 
    std::vector<TGraph*> theta_graphs; 
    std::vector<TGraph*> d0_graphs; 
    std::vector<TGraph*> z0_graphs; 
    std::vector<TH1F*> chi2_hist;  
    std::vector<TH1F*> omega_hist;  
    std::vector<TH1F*> phi_hist;  
    std::vector<TH1F*> theta_hist;  
    std::vector<TH1F*> d0_hist;  
    std::vector<TH1F*> z0_hist;  

    chi2_hist.push_back(new TH1F("chi2_1", "chi2_1", 100, 0, 50)); 
    chi2_hist.push_back(new TH1F("chi2_2", "chi2_2", 100, 0, 50)); 
    chi2_hist.push_back(new TH1F("chi2_3", "chi2_3", 100, 0, 50)); 
    chi2_hist.push_back(new TH1F("chi2_4", "chi2_4", 100, 0, 50)); 

    omega_hist.push_back(new TH1F("curvature_1", "curvature_1", 50, 0.0003, 0.00015)); 
    omega_hist.push_back(new TH1F("curvature_2", "curvature_2", 50, 0.0003, 0.00015)); 
    omega_hist.push_back(new TH1F("curvature_3", "curvature_3", 50, 0.0003, 0.00015)); 
    omega_hist.push_back(new TH1F("curvature_4", "curvature_4", 50, 0.0003, 0.00015)); 

    phi_hist.push_back(new TH1F("phi_1", "phi_1", 40, -0.05, 0.05)); 
    phi_hist.push_back(new TH1F("phi_2", "phi_2", 40, -0.05, 0.05)); 
    phi_hist.push_back(new TH1F("phi_3", "phi_3", 40, -0.05, 0.05)); 
    phi_hist.push_back(new TH1F("phi_4", "phi_4", 40, -0.05, 0.05)); 

    theta_hist.push_back(new TH1F("theta_1", "theta_1", 40, 0.02, 0.06)); 
    theta_hist.push_back(new TH1F("theta_2", "theta_2", 40, 0.02, 0.06)); 
    theta_hist.push_back(new TH1F("theta_3", "theta_3", 40, 0.02, 0.06)); 
    theta_hist.push_back(new TH1F("theta_4", "theta_4", 40, 0.02, 0.06)); 
    
    d0_hist.push_back(new TH1F("d0_1", "d0_1", 80, -10, 10)); 
    d0_hist.push_back(new TH1F("d0_2", "d0_2", 80, -10, 10)); 
    d0_hist.push_back(new TH1F("d0_3", "d0_3", 80, -10, 10)); 
    d0_hist.push_back(new TH1F("d0_4", "d0_4", 80, -10, 10)); 

    z0_hist.push_back(new TH1F("z0_1", "z0_1", 80, -2, 2)); 
    z0_hist.push_back(new TH1F("z0_2", "z0_2", 80, -2, 2)); 
    z0_hist.push_back(new TH1F("z0_3", "z0_3", 80, -2, 2)); 
    z0_hist.push_back(new TH1F("z0_4", "z0_4", 80, -2, 2)); 


    for (list<string>::iterator files_it = files.begin(); files_it != files.end(); ++files_it) { 
        t_files.push_back(new TFile((*files_it).c_str()));        
        chi2_graphs.push_back((TGraph*) t_files.back()->Get("Track #chi^{2} vs Event"));
        omega_graphs.push_back((TGraph*) t_files.back()->Get("Track #Omega vs Event"));
        phi_graphs.push_back((TGraph*) t_files.back()->Get("Track sin(#phi_{0}) vs Event"));
        theta_graphs.push_back((TGraph*) t_files.back()->Get("Track cos(#theta) vs Event"));
        d0_graphs.push_back((TGraph*) t_files.back()->Get("Track d0 vs Event"));
        z0_graphs.push_back((TGraph*) t_files.back()->Get("Track z0 vs Event"));
    }
    
    std::vector<double> chi2(chi2_graphs.size(), 0); 
    std::vector<double> omega(chi2_graphs.size(), 0); 
    std::vector<double> phi(chi2_graphs.size(), 0); 
    std::vector<double> theta(chi2_graphs.size(), 0); 
    std::vector<double> d0(chi2_graphs.size(), 0); 
    std::vector<double> z0(chi2_graphs.size(), 0); 
    std::vector<double> event(chi2_graphs.size(), 0); 

    TH1F* h_chi2_index = new TH1F("chi2_index", "chi2_index", 5, 0, 5);
    TH2F* h_z0_v_d0 = new TH2F("z0_v_d0", "z0_v_d0", 80, -2, 2, 80, -2, 2);
    TH2F* h_z0_v_chi2 = new TH2F("z0_v_chi2", "z0_v_chi2", 80, -2, 2, 25, 0, 25); 

    for (int index = 0; index < chi2_graphs.back()->GetN(); ++index) { 
        
        for (int graph_n = 0; graph_n < chi2_graphs.size(); ++graph_n) { 
            chi2_graphs[graph_n]->GetPoint(index, event[graph_n], chi2[graph_n]);
            omega_graphs[graph_n]->GetPoint(index, event[graph_n], omega[graph_n]);
            phi_graphs[graph_n]->GetPoint(index, event[graph_n], phi[graph_n]);
            theta_graphs[graph_n]->GetPoint(index, event[graph_n], theta[graph_n]);
            d0_graphs[graph_n]->GetPoint(index, event[graph_n], d0[graph_n]);
            z0_graphs[graph_n]->GetPoint(index, event[graph_n], z0[graph_n]);
        }

        double best_chi2 = 9999; 
        double best_chi2_index = -9999;
        double track_omega = 0;  

        for (int graph_n = 0; graph_n < chi2_graphs.size(); ++graph_n) {
            if (chi2[graph_n] != 0 && chi2[graph_n] < best_chi2) { 
                best_chi2 = chi2[graph_n]; 
                best_chi2_index = graph_n; 
                track_omega = omega[graph_n];
            }
        }
        
        if (best_chi2 != 9999 && best_chi2 != 0) { 
            h_chi2_index->Fill(best_chi2_index); 
            omega_hist[best_chi2_index]->Fill(omega[best_chi2_index]); 
            phi_hist[best_chi2_index]->Fill(phi[best_chi2_index]); 
            chi2_hist[best_chi2_index]->Fill(chi2[best_chi2_index]); 
            theta_hist[best_chi2_index]->Fill(theta[best_chi2_index]); 
            d0_hist[best_chi2_index]->Fill(d0[best_chi2_index]); 
            z0_hist[best_chi2_index]->Fill(z0[best_chi2_index]);
            h_z0_v_d0->Fill(z0[best_chi2_index], d0[best_chi2_index]);  
            h_z0_v_chi2->Fill(z0[best_chi2_index], chi2[best_chi2_index]);  
        }
    }

    TFile* root_file = new TFile("strategy_analysis.root", "RECREATE");

    h_chi2_index->Write();

    chi2_hist[0]->Write();  
    chi2_hist[1]->Write();  
    chi2_hist[2]->Write();  
    chi2_hist[3]->Write();

    omega_hist[0]->Write();  
    omega_hist[1]->Write();  
    omega_hist[2]->Write();  
    omega_hist[3]->Write();

    phi_hist[0]->Write();  
    phi_hist[1]->Write();  
    phi_hist[2]->Write();  
    phi_hist[3]->Write();

    theta_hist[0]->Write();  
    theta_hist[1]->Write();  
    theta_hist[2]->Write();  
    theta_hist[3]->Write();

    d0_hist[0]->Write();  
    d0_hist[1]->Write();  
    d0_hist[2]->Write();  
    d0_hist[3]->Write();

    z0_hist[0]->Write();  
    z0_hist[1]->Write();  
    z0_hist[2]->Write();  
    z0_hist[3]->Write();

    h_z0_v_d0->Write(); 
    h_z0_v_chi2->Write(); 

    delete root_file;  
}
