
#ifndef __TIMING_ANALYSIS_H__
#define __TIMING_ANALYSIS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <string>

//------------//
//--- ROOT ---//
//------------//
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>

//--------------------//
//--- HPS Analysis ---//
//--------------------//
#include <HpsAnalysis.h>
#include <Plotter.h>

//---------------//
//--- HPS DST ---//
//---------------//
#include <SvtTrack.h>
#include <SvtHit.h>

class TimingAnalysis : public HpsAnalysis { 

    public:

        /**
         * Constructor
         */
        TimingAnalysis();

        /**
         * Destructor
         */
        ~TimingAnalysis();

        /**
         *  Method used to initialize an HPS analysis.
         */
        void initialize();

        /**
         *  Method containing the code used to process an HpsEvent.
         *
         *  @param event : HpsEvent that will be processed
         */
        void processEvent(HpsEvent* event);

        /**
         *  Method used to finalize an HPS analysis.
         */
        void finalize();

        /**
         *  Method used to initialize any histograms used by the analysis.
         */
        // TODO:  This should use a histogram factory instead
        void bookHistograms();

        /**
         *  Provide a string representation of this analysis.
         *
         *  @return String representation of this analysis.
         */
        std::string toString(); 

    private:

        void fitToGaussian(TH1* histogram, double &mean, double &rms, 
                double &mean_error, double &rms_error);

        SvtTrack* track;

        Plotter* track_plotter; 
        Plotter* electron_plotter;
        Plotter* positron_plotter; 
        Plotter* top_plotter;
        Plotter* bottom_plotter;

        // Name of the class
		std::string class_name;
};

#endif
