
#ifndef __ECAL_UTILS_H__
#define __ECAL_UTILS_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <vector>

//------------------//
//   HPS Analysis   //
//------------------//
#include <FlatTupleMaker.h>

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

        /**
         *
         */
        bool hasGoodClusterPair(HpsParticle* particle); 
   
        /** */
        void setCoincidenceTime(double coin_time) { coin_time_ = coin_time; } 

        /** Use loose Ecal cluster selection. */
        void useLooseSelection(bool loose_selection) { loose_selection_ = loose_selection; } 

    private:
  
        double coin_time_{3.0}; 

        double top_time_window_low_{38.625};
        double top_time_window_high_{46.8125}; 
        double bot_time_window_low_{37.25};
        double bot_time_window_high_{47.23};

        bool loose_selection_{false}; 
};

#endif // __ECAL_UTILS_H__
