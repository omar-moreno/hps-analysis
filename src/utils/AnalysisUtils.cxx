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

double AnalysisUtils::getInvariantMass(SvtTrack* track_0, SvtTrack* track_1) {

    double p0 = AnalysisUtils::getMagnitude(track_0->getMomentum());
    double p1 = AnalysisUtils::getMagnitude(track_1->getMomentum()); 

    // Calculate the invariant mass
    double energy[2];
    double electron_mass = 0.000510998928;

    energy[0] = sqrt(p0*p0 + electron_mass*electron_mass);
    energy[1] = sqrt(p1*p1 + electron_mass*electron_mass);

    double px_sum = track_0->getMomentum()[0] + track_1->getMomentum()[0];
    double py_sum = track_0->getMomentum()[1] + track_1->getMomentum()[1];
    double pz_sum = track_0->getMomentum()[2] + track_1->getMomentum()[2];

    double p_sum = sqrt(px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);

    return sqrt(pow(energy[0]+energy[1], 2) - pow(p_sum, 2));
}

std::vector<EcalCluster*> AnalysisUtils::getClusterPair(HpsEvent* event) { 

    EcalCluster* first_cluster = NULL;
    EcalCluster* second_cluster = NULL;

    // Loop through all clusters in an event
    for (int first_cluster_n = 0; first_cluster_n < event->getNumberOfEcalClusters(); ++first_cluster_n) { 
    
        // Get an Ecal cluster from the event
        EcalCluster* cur_first_cluster = event->getEcalCluster(first_cluster_n);

        // Make sure that the Ecal cluster has a reasonable time associated 
        // with it.  If not, move on to the next cluster.
        if (cur_first_cluster->getClusterTime() <= 20) continue;

        //if (cur_first_cluster->getClusterTime() > 49 || cur_first_cluster->getClusterTime() < 40) continue;

        // Loop through the rest of the clusters and make pairs
        double min_delta_cluster_time = 1000;
        for (int second_cluster_n = (first_cluster_n + 1); second_cluster_n < event->getNumberOfEcalClusters();
                ++second_cluster_n) { 
            
            // Get another Ecal cluster from the event
            EcalCluster* cur_second_cluster = event->getEcalCluster(second_cluster_n);

            // Check if the two clusters can be considered a 'good pair'. This is done by requiring 
            // the clusters to satisfy a series of cuts
            
            // If the difference between the cluster time is greater than 2.5 ns, move on to the next cluster 
            double delta_cluster_time = 
                cur_first_cluster->getClusterTime() - cur_second_cluster->getClusterTime();
            if (std::abs(delta_cluster_time) > 2.5) continue;

            // If the difference in cluster time between the two clusters is the minimum dt in
            // the event, keep the clusters
            if (delta_cluster_time < min_delta_cluster_time) { 
                min_delta_cluster_time = delta_cluster_time;
                first_cluster = cur_first_cluster;
                second_cluster = cur_second_cluster;
            } 
        }
    }

    std::vector<EcalCluster*> pair;

    // If two pair clusters were found in the event, add them to the list of clusters.  Otherwise, 
    // return an empty pair.
    if (first_cluster != NULL && second_cluster != NULL) { 
        pair.push_back(first_cluster);
        pair.push_back(second_cluster);
    }

    return pair; 
}
