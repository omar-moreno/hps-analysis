
#include <EcalUtils.h>

EcalUtils::EcalUtils() 
    : tuple(new FlatTupleMaker("ecal_cluster_tuple.root", "results")),
      event_count(0),  
      delta_t_lower_bound(-2.0),
      delta_t_upper_bound(2.0), 
      coin_time(3.0) {
   
    this->init();
}

EcalUtils::~EcalUtils() { 
}

bool EcalUtils::hasGoodClusterPair(HpsParticle* particle) { 
   
    // Get the daughter particles composing this particle. 
    TRefArray* daughter_particles = particle->getParticles();

    // Check that the mother particle has exactly two daughters. If not, return
    // false.
    if (daughter_particles->GetEntriesFast() != 2) return false;
    
    // Check that the two daughters have an Ecal cluster associated with them.
    // If not, return false.
    if (particle->getClusters()->GetEntriesFast() != 2) return false; 
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Both daughters have clusters associated with them." << std::endl;
    //=== DEBUG
    
    std::vector<EcalCluster*> clusters = { 
        (EcalCluster*) ((HpsParticle*) particle->getClusters()->At(0)),
        (EcalCluster*) ((HpsParticle*) particle->getClusters()->At(1))
    };
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Cluster position: [" << clusters[0]->getPosition()[0] << ", " 
    //          << clusters[0]->getPosition()[1] << ", " 
    //          << clusters[0]->getPosition()[2] << "]" << std::endl; 
    //std::cout << "[ EcalUtils ]: Cluster position: [" << clusters[1]->getPosition()[0] << ", " 
    //          << clusters[1]->getPosition()[1] << ", " 
    //          << clusters[1]->getPosition()[2] << "]" << std::endl; 
    //=== DEBUG
   
    // Make sure the clusters are in opposite Ecal volumes. 
    if (clusters[0]->getPosition()[1]*clusters[1]->getPosition()[1] > 0) {
        //=== DEBUG
        //std::cout << "[ EcalUtils ]: Clusters are in same detector volume." << std::endl;
        //=== DEBUG
        return false;
    }
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Clusters are in opposite volumes." << std::endl;
    //=== DEBUG

    // Calculate the coincidence time between the two clusters.
    //double delta_cluster_time = clusters[0]->getClusterTime() - clusters[1]->getClusterTime();
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Cluster dt: " << delta_cluster_time << std::endl;
    //=== DEBUG
   
    // Check that the coincidence time between the two clusters is within some 
    // specified time window. 
    //if (abs(delta_cluster_time) >= coin_time) return false;
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Clusters are coincident" << std::endl;
    //=== DEBUG

    return true; 
}

std::vector<EcalCluster*> EcalUtils::getClusterPair(HpsEvent* event) { 

    // Initialize a vector that will be used to contain a cluster pair.  For
    // now, the clusters are set to null. 
    std::vector<EcalCluster*> cluster_pair(2, nullptr);
 
    std::vector<std::vector<EcalCluster*>> cluster_pairs; 

    // Collection used to store all clusters in an event
    std::vector<EcalCluster*> clusters; 

    // If the event has less than two clusters, return a cluster pair 
    // consisting of null pointers. 
    if (event->getNumberOfEcalClusters() < 2) return {nullptr, nullptr};
    
    // Loop over all of the clusters in the event and filter out clusters that
    // fall outside of a loose time window.  Clusters outside of this time 
    // window are very likely not to have a track associated with them.  If the
    // clusters pass the time window cut, add them to a seperate collection of
    // clusters. 
    for (int cluster_n = 0; cluster_n < event->getNumberOfEcalClusters(); ++cluster_n) { 
        double cluster_time = event->getEcalCluster(cluster_n)->getClusterTime();
        if (cluster_time < 30. || cluster_time > 60) continue;
        clusters.push_back(event->getEcalCluster(cluster_n));
    }

    // Loop through the clusters and create pairs.  Two clusters are considered
    // a pair if: 
    // 1) They are in opposite Ecal volumes i.e. top-bottom
    // 2) They are in coincidence within 1.6 ns.
    double min_delta_t = 100;
    for (int first_index = 0; first_index < clusters.size(); ++first_index) {
        
        double cluster0_y = clusters[first_index]->getPosition()[1];
        double cluster0_time = clusters[first_index]->getClusterTime();
         
        for (int second_index = (first_index + 1); second_index < clusters.size(); ++second_index) { 
       
            double cluster1_y = clusters[second_index]->getPosition()[1];
            if (cluster0_y*cluster1_y >= 0) continue; 
          
            double cluster1_time = clusters[second_index]->getClusterTime();
            double delta_t = cluster0_time - cluster1_time;
            if (abs(delta_t) > coin_time) continue;
            
            if (delta_t > min_delta_t) continue;
           
            min_delta_t = delta_t;
            cluster_pair[0] = clusters[first_index];
            cluster_pair[1] = clusters[second_index]; 
        } 
    }
    
    //tuple->fill(); 

    return cluster_pair;     
}

void EcalUtils::init() { 

    tuple->addVariable("cluster0_energy");
    tuple->addVariable("cluster1_energy");
    tuple->addVariable("cluster0_x");
    tuple->addVariable("cluster0_y");
    tuple->addVariable("cluster1_x");
    tuple->addVariable("cluster1_y");
    tuple->addVariable("cluster0_time");
    tuple->addVariable("cluster1_time");
    tuple->addVariable("event_n"); 
    tuple->addVariable("n_clusters");
    tuple->addVariable("n_pairs"); 
    tuple->addVariable("pass_cluster_n_cut");

    tuple->addVector("cluster_time");
    tuple->addVector("cluster_pair_dt");
    tuple->addVector("cluster_pair_energy_diff"); 
    tuple->addVector("cluster_pair_energy_sum");
    tuple->addVector("pass_cluster_time_cut");
    tuple->addVector("pass_geo_cut");
    tuple->addVector("pass_coin_cut");
}

/*
void EcalUtils::bookHistograms() { 
  
    //   Ecal cluster pair energy   //
    //------------------------------//

    plotter->build2DHistogram("cluster pair energy", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy")->GetXaxis()->SetTitle("Ecal cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy")->GetYaxis()->SetTitle("Ecal cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetXaxis()->SetTitle("Ecal cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->GetYaxis()->SetTitle("Ecal cluster energy (GeV)");

    plotter->build2DHistogram("cluster pair energy - cuts: fiducial, time", 50, 0, 1.5, 50, 0, 1.5);
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, time")->GetXaxis()->SetTitle("Ecal cluster energy (GeV)");
    plotter->get2DHistogram("cluster pair energy - cuts: fiducial, time")->GetYaxis()->SetTitle("Ecal cluster energy (GeV)");

    plotter->build1DHistogram("cluster energy sum", 100, 0, 1.5)->GetXaxis()->SetTitle("Ecal cluster energy sum (GeV)"); 
    plotter->build1DHistogram("cluster energy sum - cuts: fiducial", 100, 0, 1.5)->GetXaxis()->SetTitle("Ecal cluster energy sum (GeV)");
    plotter->build1DHistogram("cluster energy sum - cuts: fiducial, time", 100, 0, 1.5)->GetXaxis()->SetTitle("Ecal cluster energy sum (GeV)"); 
    
    //   Ecal cluster pair time   //
    //----------------------------//

    plotter->build1DHistogram("cluster pair #Delta t", 80, -10, 10)->GetXaxis()->SetTitle("Cluster pair #Delta t (ns)");
    plotter->build1DHistogram("cluster pair #Delta t - cuts: fiducial", 80, -10, 10)->GetXaxis()->SetTitle("Cluster pair #Delta t (ns)");
    plotter->build1DHistogram("cluster pair #Delta t - cuts: fiducial, time", 80, -10, 10)->GetXaxis()->SetTitle("Cluster pair #Delta t (ns)");

    plotter->build2DHistogram("cluster pair time", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time")->GetXaxis()->SetTitle("Ecal cluster time (ns)");
    plotter->get2DHistogram("cluster pair time")->GetYaxis()->SetTitle("Ecal cluster time (ns)");

    plotter->build2DHistogram("cluster pair time - cuts: fiducial", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetXaxis()->SetTitle("Ecal cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial")->GetYaxis()->SetTitle("Ecal cluster time (ns)");
    
    plotter->build2DHistogram("cluster pair time - cuts: fiducial, time", 250, 0, 125, 250, 0, 125);
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, time")->GetXaxis()->SetTitle("Ecal cluster time (ns)");
    plotter->get2DHistogram("cluster pair time - cuts: fiducial, time")->GetYaxis()->SetTitle("Ecal cluster time (ns)");

    //   Ecal cluster position   //
    //---------------------------//

    plotter->build1DHistogram("cluster pair delta x", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair delta x (mm)");
    plotter->build1DHistogram("cluster pair delta x - cuts: fiducial", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair delta x (mm)");
    plotter->build1DHistogram("cluster pair delta x - cuts: fiducial, time", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair delta x (mm)");

    plotter->build1DHistogram("cluster pair x sum", 200, -300, -100)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
    plotter->build1DHistogram("cluster pair x sum - cuts: fiducial", 200, -300, -100)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
    plotter->build1DHistogram("cluster pair x sum - cuts: fiducial, time", 200, -300, 100)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");

    plotter->build2DHistogram("cluster x vs cluster x", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y vs cluster y", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y vs cluster y")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y vs cluster y")->GetYaxis()->SetTitle("Cluster position - y (mm)");

    plotter->build2DHistogram("cluster x vs cluster x - cuts: fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y vs cluster y - cuts: fiducial", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    
    plotter->build2DHistogram("cluster x vs cluster x - cuts: fiducial, time", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, time")->GetYaxis()->SetTitle("Cluster position - x (mm)");
    plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, time")->GetYaxis()->SetTitle("Cluster position - x (mm)");

    plotter->build2DHistogram("cluster y vs cluster y - cuts: fiducial, time", 200, -200, 200, 200, -200, 200);
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, time")->GetYaxis()->SetTitle("Cluster position - y (mm)");
    plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, time")->GetYaxis()->SetTitle("Cluster position - y (mm)");

}*/
