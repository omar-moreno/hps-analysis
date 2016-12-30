
#include <EcalUtils.h>

EcalUtils::EcalUtils() { 
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

    double top_index = 0;
    double bot_index = 1;
    if (((EcalCluster*) ((HpsParticle*) particle->getClusters()->At(bot_index)))->getPosition()[1] > 0) { 
        top_index = 1;
        bot_index = 0;
    }

    EcalCluster* top{(EcalCluster*) ((HpsParticle*) particle->getClusters()->At(top_index))};
    EcalCluster* bot{(EcalCluster*) ((HpsParticle*) particle->getClusters()->At(bot_index))};

    // Make sure the clusters are in opposite Ecal volumes. 
    if (top->getPosition()[1]*bot->getPosition()[1] > 0) {
        return false;
    }

    if (loose_selection_) return true;

    if (top->getClusterTime() < top_time_window_low_ || top->getClusterTime() > top_time_window_high_) return false; 
    if (bot->getClusterTime() < bot_time_window_low_ || bot->getClusterTime() > bot_time_window_high_) return false; 
    
    double coin_time = top->getClusterTime() - bot->getClusterTime();
    if (abs(coin_time) > coin_time_) return false;

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
    /*double min_delta_t = 100;
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
    }*/
    
    //tuple->fill(); 

    return cluster_pair;     
}

