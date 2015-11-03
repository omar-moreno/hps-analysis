'''
'''

def is_good_data_event(event) : 

    if not event.isSingle1Trigger() : return False

    if not event.isSvtBiasOn() : return False

    if not event.isSvtClosed() : return False

    if event.hasSvtBurstModeNoise() : return False

    if event.hasSvtEventHeaderErrors() : return False

    return True

def get_good_cluster_pair(event) :
 
    cluster_pair = []

    for first_cluster_n in xrange(event.getNumberOfEcalClusters()) :
        
        first_cluster = event.getEcalCluster(first_cluster_n)

        for second_cluster_n in xrange(event.getNumberOfEcalClusters()) :

            second_cluster = event.getEcalCluster(second_cluster_n)

            if first_cluster.getPosition()[1]*second_cluster.getPosition()[1] > 0 : continue

            cluster_pair_dt = first_cluster.getClusterTime() - second_cluster.getClusterTime()
                
            if cluster_pair_dt < -1.6 or cluster_pair_dt > 1.7 : continue

            cluster_pair = [first_cluster, second_cluster]

    return cluster_pair
