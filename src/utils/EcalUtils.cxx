
#include <EcalUtils.h>

EcalUtils::EcalUtils() 
    : plotter(new Plotter()), 
      delta_t_lower_bound(-1.6),
      delta_t_upper_bound(1.7) {
    
    this->bookHistograms(); 
}

EcalUtils::~EcalUtils() { 
    delete plotter; 
}

bool EcalUtils::isGoodClusterPair(HpsParticle* particle) { 
    
    TRefArray* cluster_objects = particle->getClusters();
    if (cluster_objects->GetEntriesFast() != 2) return false;

    std::vector<EcalCluster*> clusters = { 
        (EcalCluster*) cluster_objects->At(0), 
        (EcalCluster*) cluster_objects->At(1)
    }; 

    if (clusters[0]->getPosition()[1]*clusters[1]->getPosition()[1] > 0) return false;

    double delta_cluster_time 
        = clusters[0]->getClusterTime() - clusters[1]->getClusterTime();

    if (delta_cluster_time < delta_t_lower_bound || delta_cluster_time > delta_t_upper_bound) return false;

    return true; 
}

std::vector<EcalCluster*> EcalUtils::getClusterPair(HpsEvent* event) { 

    std::vector<EcalCluster*> cluster_pair(2, nullptr); 

    // Loop through all clusters in an event and find the best pair
    for (int first_cluster_n = 0; first_cluster_n < event->getNumberOfEcalClusters(); ++first_cluster_n) { 
    
        // Get an Ecal cluster from the event
        EcalCluster* current_first_cluster = event->getEcalCluster(first_cluster_n);

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
            if (current_first_cluster->getPosition()[1]*current_second_cluster->getPosition()[1] > 0) continue;

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

    plotter->build1DHistogram("cluster pair x sum", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
    plotter->build1DHistogram("cluster pair x sum - cuts: fiducial", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");
    plotter->build1DHistogram("cluster pair x sum - cuts: fiducial, time", 200, -200, 200)->GetXaxis()->SetTitle("Ecal cluster pair x sum (mm)");

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
