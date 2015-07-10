
#include <TimingAnalysis.h>

TimingAnalysis::TimingAnalysis()
    : track(NULL),
      track_plotter(new Plotter()),
      electron_plotter(new Plotter()),
      positron_plotter(new Plotter()),
      top_plotter(new Plotter()), 
      bottom_plotter(new Plotter()), 
      class_name("TimingAnalysis") {  
}

TimingAnalysis::~TimingAnalysis() { 
    delete track_plotter; 
    delete electron_plotter;
    delete positron_plotter; 
    delete top_plotter;
    delete bottom_plotter;
}

void TimingAnalysis::initialize() { 
    this->bookHistograms();
}

void TimingAnalysis::processEvent(HpsEvent* event) {

    //if (!event->isPair1Trigger()) return;
    if (!event->isSingle1Trigger()) return;

    // Loop over all of the tracks in the event
    for (int track_n = 0; track_n < event->getNumberOfTracks(); ++track_n) { 

        // Get a track from the event
        track = event->getTrack(track_n); 

        if (track->isTopTrack()) {
            track_plotter->get1DHistogram("Top track time")->Fill(track->getTrackTime());
        } else { 
            track_plotter->get1DHistogram("Bottom track time")->Fill(track->getTrackTime());
        }
    
    
        TRefArray* hits = track->getSvtHits();
        double total_hits = (double) hits->GetSize();
        for (int hit_n = 0; hit_n < hits->GetSize(); ++hit_n) { 

            int layer = ((SvtHit*) hits->At(hit_n))->getLayer();
            double t0_residual = track->getTrackTime() - ((SvtHit*) hits->At(hit_n))->getTime();
            t0_residual = t0_residual/sqrt((total_hits - 1)/total_hits);
            
            if (track->isTopTrack()) {
                std::string top_name = "Top Layer " + std::to_string(layer);
                track_plotter->get1DHistogram(top_name + " - Track time - Hit time")->Fill(t0_residual);
            } else { 
                std::string bot_name = "Bottom Layer " + std::to_string(layer);
                track_plotter->get1DHistogram(bot_name + " - Track time - Hit time")->Fill(t0_residual);
            } 
        } 
    }
}

void TimingAnalysis::finalize() {

    TGraphErrors* top_t0 = (TGraphErrors*) track_plotter->buildGraph("SVT top t0 resolution");
    TGraphErrors* bot_t0 = (TGraphErrors*) track_plotter->buildGraph("SVT bottom t0 resolution");

    double mean = 0;
    double mean_error = 0; 
    double rms = 0;
    double rms_error = 0;  
    for (int module_n = 1; module_n <= 6; ++module_n) {
        std::string top_name = "Top Layer " + std::to_string(module_n);
        std::string bot_name = "Bottom Layer " + std::to_string(module_n);
        
        this->fitToGaussian(track_plotter->get1DHistogram(top_name + " - Track time - Hit time"), mean, rms, mean_error, rms_error);
        top_t0->SetPoint(module_n - 1, module_n, rms);
        top_t0->SetPointError(module_n - 1, module_n, rms_error);
        

        this->fitToGaussian(track_plotter->get1DHistogram(bot_name + " - Track time - Hit time"), mean, rms, mean_error, rms_error);
        bot_t0->SetPoint(module_n - 1, module_n, rms);
        bot_t0->SetPointError(module_n - 1, module_n, rms_error);
    }

    track_plotter->saveToRootFile("timing_analysis.root");
    track_plotter->saveToPdf("timing_analysis.pdf");
}

void TimingAnalysis::bookHistograms() { 

    track_plotter->build1DHistogram("Top track time", 100, -10, 10)->GetXaxis()->SetTitle("Track time [ns]");
    track_plotter->build1DHistogram("Bottom track time", 100, -10, 10)->GetXaxis()->SetTitle("Track time [ns]");


    for (int module_n = 1; module_n <= 6; ++module_n) {
        std::string top_name = "Top Layer " + std::to_string(module_n);
        std::string bot_name = "Bottom Layer " + std::to_string(module_n);
    
        track_plotter->build1DHistogram(top_name + " - Track time - Hit time", 100, -10, 10)->GetXaxis()->SetTitle("Track time - Hit time [ns]");
        track_plotter->build1DHistogram(bot_name + " - Track time - Hit time", 100, -10, 10)->GetXaxis()->SetTitle("Track time - Hit time [ns]");
    
    }
}


void TimingAnalysis::fitToGaussian(TH1* histogram, double &mean, double &rms, 
       double &rms_error, double &mean_error ) { 
    
    mean = histogram->GetMean(); 
    rms = histogram->GetRMS();

    TF1* gaussian = new TF1("gaussian", "gaus");
    gaussian->SetRange(mean - 3*rms, mean + 3*rms);
    histogram->Fit("gaussian", "RQ");
    mean = gaussian->GetParameter(1);
    mean_error = gaussian->GetParError(1);
    rms = gaussian->GetParameter(2);
    rms_error = gaussian->GetParError(2);
    gaussian->SetRange(mean - 5*rms, mean + 5*rms);
    delete gaussian; 
}

std::string TimingAnalysis::toString() { 
    std::string string_buffer = "Class Name: " + class_name; 
    return string_buffer;   
}

