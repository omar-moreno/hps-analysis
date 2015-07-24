/**
 * @file AnalysisUtils.cxx
 * @brief A set of utilities commonly used when doing analysis
 * @author Omar Moreno <omoreno1@ucsc.edu> \n
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date May 20, 2015
 */

#include <AnalysisUtils.h>

double AnalysisUtils::getMagnitude(std::vector<double> v) { 
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

std::vector<EcalCluster*> AnalysisUtils::getClusterPair(HpsEvent* event) { 

    EcalCluster* first_cluster = NULL;
    EcalCluster* second_cluster = NULL;

    // Loop through all clusters in an event
    for (int first_cluster_n = 0; first_cluster_n < event->getNumberOfEcalClusters(); ++first_cluster_n) { 
    
        // Get an Ecal cluster from the event
        EcalCluster* cur_first_cluster = event->getEcalCluster(first_cluster_n);

        if (cur_first_cluster->getClusterTime() > 49 || cur_first_cluster->getClusterTime() < 40) continue;

        // Loop through the rest of the clusters and make pairs
        for (int second_cluster_n = (first_cluster_n + 1); second_cluster_n < event->getNumberOfEcalClusters();
                ++second_cluster_n) { 
            
            // Get another Ecal cluster from the event
            EcalCluster* cur_second_cluster = event->getEcalCluster(second_cluster_n);

            // Check if the two clusters can be considered a 'good pair'. This is done by requiring 
            // the clusters to satisfy a series of cuts
             
            double delta_cluster_time = 
                cur_first_cluster->getClusterTime() - cur_second_cluster->getClusterTime();
            if (std::abs(delta_cluster_time) > 2.5) continue;

            if (cur_second_cluster->getClusterTime() > 50 || cur_second_cluster->getClusterTime() < 39) continue;

            first_cluster = cur_first_cluster;
            second_cluster = cur_second_cluster;
            break;
        }
    }

    std::vector<EcalCluster*> pair;

    if (first_cluster != NULL && second_cluster != NULL) { 
        pair.push_back(first_cluster);
        pair.push_back(second_cluster);
    }

    return pair; 
}
