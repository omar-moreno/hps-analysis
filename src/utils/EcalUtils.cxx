
#include <EcalUtils.h>

EcalUtils::EcalUtils() 
    : plotter(new Plotter()), 
      delta_t_lower_bound(-1.6),
      delta_t_upper_bound(1.6) {
    
    this->bookHistograms(); 
}

EcalUtils::~EcalUtils() { 
    delete plotter; 
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
    
    if (clusters[0]->getPosition()[1]*clusters[1]->getPosition()[1] > 0) {
        //=== DEBUG
        //std::cout << "[ EcalUtils ]: Clusters are in same detector volume." << std::endl;
        //=== DEBUG
        return false;
    }
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Clusters are in opposite volumes." << std::endl;
    //=== DEBUG

    double delta_cluster_time = clusters[0]->getClusterTime() - clusters[1]->getClusterTime();
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Cluster dt: " << delta_cluster_time << std::endl;
    //=== DEBUG
    
    if (delta_cluster_time < delta_t_lower_bound || delta_cluster_time > delta_t_upper_bound) return false;
    
    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Clusters are coincident" << std::endl;
    //=== DEBUG

    return true; 
}

std::vector<EcalCluster*> EcalUtils::getClusterPair(HpsEvent* event) { 

    // Initialize a vector that will be used to contain a cluster pair.  For
    // now, the clusters are set to null. 
    std::vector<EcalCluster*> cluster_pair(2, nullptr); 
    
    // If the event has less than two clusters, return a cluster pair 
    // consisting of null pointers. 
    if (event->getNumberOfEcalClusters() < 2) return cluster_pair; 

    //=== DEBUG
    //std::cout << "[ EcalUtils ]: Searching for best cluster pair " << std::endl;
    //std::cout << "[ EcalUtils ]: Number of clusters: " << event->getNumberOfEcalClusters() << std::endl;
    //=== DEBUG

    // Loop through all clusters in an event and find the best pair
    for (int first_cluster_n = 0; first_cluster_n < event->getNumberOfEcalClusters(); ++first_cluster_n) { 
    
        // Get an Ecal cluster from the event
        EcalCluster* current_first_cluster = event->getEcalCluster(first_cluster_n);

        //if (first_cluster_n == 0) { 
            //std::cout << "[ EcalUtils ]: First cluster number: " << first_cluster_n << std::endl;
            //std::cout << "[ EcalUtils ]: First cluster time: " << current_first_cluster->getClusterTime() << std::endl;
            //    std::cout << "[ EcalUtils ]: Ecal cluster position: [ " 
            //        << current_first_cluster->getPosition()[0] << ", " 
            //        << current_first_cluster->getPosition()[1] << ", " 
             //       << current_first_cluster->getPosition()[2] << " ]" 
            //        << std::endl;
        //}

        // Make sure that the Ecal cluster has a reasonable time associated 
        // with it.  If not, move on to the next cluster.
        //if (cur_first_cluster->getClusterTime() <= 20) continue;

        // Loop through the rest of the clusters in the event and check if 
        // a match can be found
        double min_delta_cluster_time = 1000;
        for (int second_cluster_n = (first_cluster_n + 1); second_cluster_n < event->getNumberOfEcalClusters();
                ++second_cluster_n) { 


            // Get the second Ecl cluster in the event
            EcalCluster* current_second_cluster = event->getEcalCluster(second_cluster_n);
        
            //if (first_cluster_n == 0) { 
                //std::cout << "[ EcalUtils ]: Second cluster number: " << second_cluster_n << std::endl;
                //std::cout << "[ EcalUtils ]: Second cluster time: " << current_second_cluster->getClusterTime() << std::endl;
                //std::cout << "[ EcalUtils ]: Ecal cluster position: [ " 
                //    << current_second_cluster->getPosition()[0] << ", " 
                //    << current_second_cluster->getPosition()[1] << ", " 
                //    << current_second_cluster->getPosition()[2] << " ]" 
                //    << std::endl;
            
            //}

            plotter->get2DHistogram("cluster pair energy")->Fill(current_first_cluster->getEnergy(), 
                    current_second_cluster->getEnergy());
            plotter->get1DHistogram("cluster energy sum")->Fill(current_first_cluster->getEnergy()
                   + current_second_cluster->getEnergy());  

            // If the difference between the cluster time is greater than 2.5 ns, move on to the next cluster 
            double delta_cluster_time 
                = current_first_cluster->getClusterTime() - current_second_cluster->getClusterTime();

            plotter->get1DHistogram("cluster pair #Delta t")->Fill(delta_cluster_time);
            plotter->get2DHistogram("cluster pair time")->Fill(current_first_cluster->getClusterTime(), 
                    current_second_cluster->getClusterTime()); 
            plotter->get2DHistogram("cluster x vs cluster x")->Fill(current_first_cluster->getPosition()[0],
                    current_second_cluster->getPosition()[0]);
            plotter->get2DHistogram("cluster y vs cluster y")->Fill(current_first_cluster->getPosition()[1],
                    current_second_cluster->getPosition()[1]);
            plotter->get1DHistogram("cluster pair delta x")->Fill(current_first_cluster->getPosition()[0] 
                    - current_second_cluster->getPosition()[0]); 
            plotter->get1DHistogram("cluster pair x sum")->Fill(current_first_cluster->getPosition()[0] 
                    + current_second_cluster->getPosition()[0]); 

            // Require that the two clusters are in opposite volumes.  This cut should
            // eventually become part of the standard pair requirement.
            if (current_first_cluster->getPosition()[1]*current_second_cluster->getPosition()[1] > 0) { 
                if (first_cluster_n == 0) { 
                    //std::cout << "[ EcalUtils ]: Clusters are in same volume." << std::endl;
                }
                continue;
            }

            plotter->get2DHistogram("cluster pair energy - cuts: fiducial")->Fill(current_first_cluster->getEnergy(), 
                    current_second_cluster->getEnergy()); 
            plotter->get1DHistogram("cluster pair #Delta t - cuts: fiducial")->Fill(delta_cluster_time);
            plotter->get2DHistogram("cluster pair time - cuts: fiducial")->Fill(current_first_cluster->getClusterTime(), 
                    current_second_cluster->getClusterTime()); 
            plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial")->Fill(current_first_cluster->getPosition()[0],
                    current_second_cluster->getPosition()[0]);
            plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial")->Fill(current_first_cluster->getPosition()[1],
                    current_second_cluster->getPosition()[1]);
            plotter->get1DHistogram("cluster energy sum - cuts: fiducial")->Fill(current_first_cluster->getEnergy()
                   + current_second_cluster->getEnergy());  
            plotter->get1DHistogram("cluster pair delta x - cuts: fiducial")->Fill(current_first_cluster->getPosition()[0] 
                    - current_second_cluster->getPosition()[0]); 
            plotter->get1DHistogram("cluster pair x sum - cuts: fiducial")->Fill(current_first_cluster->getPosition()[0] 
                    + current_second_cluster->getPosition()[0]); 

            if (delta_cluster_time < delta_t_lower_bound || delta_cluster_time > delta_t_upper_bound) continue;

            plotter->get2DHistogram("cluster pair energy - cuts: fiducial, time")->Fill(current_first_cluster->getEnergy(), 
                    current_second_cluster->getEnergy()); 
            plotter->get1DHistogram("cluster pair #Delta t - cuts: fiducial, time")->Fill(delta_cluster_time);
            plotter->get2DHistogram("cluster pair time - cuts: fiducial, time")->Fill(current_first_cluster->getClusterTime(), 
                    current_second_cluster->getClusterTime()); 
            plotter->get2DHistogram("cluster x vs cluster x - cuts: fiducial, time")->Fill(current_first_cluster->getPosition()[0],
                    current_second_cluster->getPosition()[0]);
            plotter->get2DHistogram("cluster y vs cluster y - cuts: fiducial, time")->Fill(current_first_cluster->getPosition()[1],
                    current_second_cluster->getPosition()[1]);
            plotter->get1DHistogram("cluster energy sum - cuts: fiducial, time")->Fill(current_first_cluster->getEnergy()
                   + current_second_cluster->getEnergy());  
            plotter->get1DHistogram("cluster pair delta x - cuts: fiducial, time")->Fill(current_first_cluster->getPosition()[0] 
                    - current_second_cluster->getPosition()[0]); 
            plotter->get1DHistogram("cluster pair x sum - cuts: fiducial, time")->Fill(current_first_cluster->getPosition()[0] 
                    + current_second_cluster->getPosition()[0]); 

            // If the difference in cluster time between the two clusters is the minimum dt in
            // the event, keep the clusters
            if (delta_cluster_time < min_delta_cluster_time) { 
                min_delta_cluster_time = delta_cluster_time;
                cluster_pair[0] = current_first_cluster; 
                cluster_pair[1] = current_second_cluster; 
            } 
        }
    }
    
    if (cluster_pair[0] != nullptr && cluster_pair[1] != nullptr) {
        //std::cout << "[ EcalUtils ]: cluster 0 time : " << cluster_pair[0]->getClusterTime() << " cluster 1 time: " << cluster_pair[1]->getClusterTime() << std::endl;
        //std::cout << "[ EcalUtils ]: Ecal cluster 0 position: [ " << cluster_pair[0]->getPosition()[0] << ", " << cluster_pair[0]->getPosition()[1] << ", " 
        //          << cluster_pair[0]->getPosition()[2] << " ]" << std::endl;
        //std::cout << "[ EcalUtils ]: Ecal cluster 1 position: [ " << cluster_pair[1]->getPosition()[0] << ", " << cluster_pair[1]->getPosition()[1] << ", " 
        //          << cluster_pair[1]->getPosition()[2] << " ]" << std::endl;
    }
    return cluster_pair; 
}

void EcalUtils::saveHistograms() { 
    plotter->saveToRootFile("ecal_cluster_plots.root");  
}

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

}
