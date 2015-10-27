from MollerPreprocessor import MollerPreprocessor

class MollerBkgPreprocessor(MollerPreprocessor) : 
    
    def is_good_cluster_pair(self, cluster_pair) : 

        if (super(MollerBkgPreprocessor, self).is_good_cluster_pair(cluster_pair) 
                and (cluster_pair[0].getEnergy() > 0.8 or cluster_pair[1].getEnergy() > 0.8)) : 
            return True

        return False
