
/** 
 * @file BumpHunter.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date January 14, 2015
 *
 */

#include <BumpHunter.h>

BumpHunter::BumpHunter(int poly_order) 
    : comp_model(nullptr), 
      bkg_model(nullptr),
      model(nullptr),
      signal(nullptr), 
      bkg(nullptr),
      ofs(nullptr), 
      window_size(0.01),
      bkg_poly_order(poly_order), 
      bkg_only(false) {

    // Turn off all messages except errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    // Independent variable
    variable_map["invariant mass"] = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0., 0.1);

    //   Signal PDF   //
    //----------------//   

    variable_map["A' mass"]  = new RooRealVar("A' Mass",  "A' Mass",  0.03);
    variable_map["A' mass resolution"] 
        = new RooRealVar("A' Mass Resolution", "A' Mass Resolution", this->getMassResolution(0.03));

    signal = new RooGaussian("signal", "signal", *variable_map["invariant mass"],
                             *variable_map["A' mass"], *variable_map["A' mass resolution"]);

    //   Bkg PDF   //
    //-------------//

    std::string name;
    for (int order = 1; order <= bkg_poly_order; ++order) {
        name = "t" + std::to_string(order);
        variable_map[name] = new RooRealVar(name.c_str(), name.c_str(), 0, -2, 2);
        arg_list.add(*variable_map[name]);
    } 

    bkg = new RooChebychev("bkg", "bkg", *variable_map["invariant mass"], arg_list);

    //   Composite Models   //
    //----------------------//

    variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -10000, 10000);
    variable_map["bkg yield"] = new RooRealVar("bkg yield", "bkg yield", 300000, 100, 10000000);

    comp_model = new RooAddPdf("comp model", "comp model", RooArgList(*signal, *bkg), 
                               RooArgList(*variable_map["signal yield"], *variable_map["bkg yield"]));

    bkg_model = new RooAddPdf("bkg model", "bkg model", 
                              RooArgList(*bkg), RooArgList(*variable_map["bkg yield"]));
    model = comp_model;
}


BumpHunter::~BumpHunter() {

    for (auto& element : variable_map) { 
       delete element.second; 
    }
    variable_map.clear();

    delete signal;
    delete bkg;
    delete comp_model; 
    if (ofs != nullptr) ofs->close();
}

std::map<double, RooFitResult*> BumpHunter::fit(TH1* histogram, double start, double end, double window_step) {

    // Set the range of the mass variable based on the range of the histogram. 
    variable_map["invariant mass"]->setRange(histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax()); 

    // Create a histogram object compatible with RooFit.
    RooDataHist* data = new RooDataHist("data", "data", RooArgList(*variable_map["invariant mass"]), histogram);
    
    // Create a container for the results from the fit to each window.
    std::map<double, RooFitResult*> results; 
    
    while (start <= (end - window_size)) { 
     
        // Fit the histogram within a window (start, start + window_size) and
        // save the result to the map of results.
        double ap_hypothesis = start + window_size/2;
        RooFitResult* result = this->fit(data, start);
        results[ap_hypothesis] = result; 
        if (ofs != nullptr) result->printMultiline(*ofs, 0, kTRUE, "");
        
        start += window_step; 
        
        this->resetParameters(result->floatParsInit()); 
    }
    
    delete data;

    return results;  
} 

/*
std::map<double, RooFitResult*> BumpHunter::fit(TH1* histogram, double start, double end, double window_step) { 
        
        if (index == 1 && test) { 
            RooPlot* frame = variable_map["signal yield"]->frame(RooFit::Range(-1000, 1000)); 
            //nll->plotOn(frame, RooFit::ShiftToZero());
            auto nllp = nll->createProfile(*variable_map["signal yield"]); 
            nllp->plotOn(frame, RooFit::Name("nll"));
            RooCurve* curve = (RooCurve*) frame->getCurve("nll"); 
            std::cout << "x: " << 0 << " nll: " << curve->findPoint(0) << std::endl;

            TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600); 
            frame->Draw(); 
            canvas->SaveAs("ll_profile2.C"); 
            delete canvas; 
            test = false;  
        }
}*/

RooFitResult* BumpHunter::fit(TH1* histogram, double window_start) { 

    // Set the range of the mass variable based on the range of the histogram. 
    variable_map["invariant mass"]->setRange(histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax()); 

    // Create a histogram object compatible with RooFit.
    RooDataHist* data = new RooDataHist("data", "data", RooArgList(*variable_map["invariant mass"]), histogram);

    // Fit the result within a window (window_start, window_start + window_size)
    RooFitResult* result = this->fit(data, window_start);
    
    // Delete the histogram object from memory
    delete data;

    // Return the result
    return result;
}

RooFitResult* BumpHunter::fit(RooDataHist* data, double window_start) { 

    // Set the A' mass hypothesis to the middle of the windowd 
    // Should there be a check to make sure the A' mass hypothesis is not set to
    // the edge of a bin?
    double ap_hypothesis = (window_start + window_size/2); 
    variable_map["A' mass"]->setVal(ap_hypothesis);
    variable_map["A' mass resolution"]->setVal(this->getMassResolution(ap_hypothesis)); 
    
    // Set the range that will be used in the fit
    std::string range_name = "ap_mass_" + std::to_string(ap_hypothesis); 
    variable_map["invariant mass"]->setRange(range_name.c_str(), window_start, window_start + window_size); 
   
    // Estimate the background normalization within the window by integrating
    // the histogram within the window range.  This should be close to the 
    // background yield in the case where there is no signal present.
    double integral = data->sumEntries(0, range_name.c_str()); 
    variable_map["bkg yield"]->setVal(integral);
    //if (ofs != nullptr) ofs << "Estimated bkg in range (" << start << ", " << start + window_size << "): " << integral; 
    //<< std::endl;

    // Construct a log likelihood using the data set within the range specified
    // above.  This is equivalent to saying that 
    // nll = -ln(L( window_start < x < window_start + window_size | mu, theta))
    // where mu is the signal yield and theta represents the set of all 
    // nuisance parameters which in this case are the background normalization
    // and polynomial constants.  Since the likelihood is being constructed
    // from a histogram, use an extended likelihood.
    RooAbsReal* nll = model->createNLL(*data, 
            RooFit::Extended(kTRUE), 
            RooFit::SumCoefRange(range_name.c_str()), 
            RooFit::Range(range_name.c_str()),
            RooFit::Verbose(kFALSE),
            RooFit::NumCPU(2));

    // Instantiate minuit using the constructed likelihood above
    RooMinuit m(*nll);

    // Turn off all print out
    m.setPrintLevel(-1000);


    // Use migrad to minimize the likelihood.  If migrad fails to find a minimum,
    // run simplex in order to run a sparser search for a minimum followed by
    // migrad again.
    int status = m.migrad(); 
    if (status != 0) { 
        m.simplex();
        status = m.migrad();
    }

    // If a valid minimum was found, then
    // 1) Run improve to try and find a better minimum.  This is done in case 
    //    migrad ends up finding a local minium instead of a global one.
    // 2) Run hesse
    // 3) Run minos in order to optimize the errors the signal yield.
    if (status == 0) { 
        m.improve();
        m.hesse();
        m.minos(*variable_map["signal yield"]);
    }

    // Save the results of the fit
    RooFitResult* result = m.save(); 

    // Delete the constructed negative log likelihood
    delete nll; 

    // Return the saves result
    return result; 
} 

void BumpHunter::resetParameters(RooArgList initial_params) { 
    
    for (auto& element : variable_map) { 
        if (initial_params.find(element.second->GetName()) == NULL) continue;

        RooRealVar* var = (RooRealVar*) initial_params.find(element.second->GetName());

        element.second->setVal(var->getVal());
        element.second->setError(var->getError()); 

    }
}

void BumpHunter::fitBkgOnly() { 
    this->bkg_only = true; 
    model = bkg_model; 
}

void BumpHunter::writeResults() { 
    
    // Create the output file name string
    char buffer[100];
    sprintf(buffer, "results_order%i_window%i.txt", bkg_poly_order, window_size*1000);

    // Create a file stream  
    ofs = new std::ofstream(buffer, std::ofstream::out); 
}
