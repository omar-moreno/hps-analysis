
#ifndef __ECAL_UTILS_H__
#define __ECAL_UTILS_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <vector>

//------------------//
//   HPS Analysis   //
//------------------//
#include <Plotter.h>   

//-------------//
//   HPS DST   //
//-------------//
#include <HpsEvent.h>
#include <EcalCluster.h>

class EcalUtils { 

    public: 

        /** Constructor */
        EcalUtils(); 

        /** Destructor */
        ~EcalUtils(); 

        /**
         * Loop through all Ecal clusters in an event and find a clusters that
         * form a pair.
         *
         * @param event HpsEvent object used to retrieve Ecal clusters in an 
         *              event.
         * @return A vector containing a pair of clusters.
         */
        // TODO: This should return a std::pair
        std::vector<EcalCluster*> getClusterPair(HpsEvent* event);

        /** Save the histogram to a ROOT file. */
        void saveHistograms();

    private:

        /** Book histograms */
       void bookHistograms();

        /** Plotter */
        Plotter* plotter;  

        /** */
        double delta_t_lower_bound;
        
        /** */
        double delta_t_upper_bound; 
};

#endif // __ECAL_UTILS_H__
