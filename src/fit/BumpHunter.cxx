
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
      low_bound(-9999), 
      high_bound(-9999),
      max_window_size(0.02205),
      window_size(0.01),
      bkg_poly_order(poly_order), 
      bkg_only(false),
      debug(true), 
      fix_window(false) {

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

    variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -100000, 100000);
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

std::map<double, HpsFitResult*> BumpHunter::fitWindow(TH1* histogram, double start, double end, double step) {

    // Set the range of the mass variable based on the range of the histogram. 
    variable_map["invariant mass"]->setRange(histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax()); 

    // Create a histogram object compatible with RooFit.
    RooDataHist* data = new RooDataHist("data", "data", RooArgList(*variable_map["invariant mass"]), histogram);
    
    // Create a container for the results from the fit to each window.
    std::map<double, HpsFitResult*> results; 
    
    while (start <= end) { 
     
        // Construct a window of size x*(A' mass resolution) around the A' mass
        // hypothesis and do a Poisson likelihood fit within the window range. 
        HpsFitResult* result = this->fitWindow(data, start);

        // Save the result to the map of results
        results[start] = result; 
        if (ofs != nullptr) result->getRooFitResult()->printMultiline(*ofs, 0, kTRUE, "");

        // Increment the A' mass hypothesis
        start += step; 
    
        // Reset all of the parameters to their original values
        this->resetParameters(result->getRooFitResult()->floatParsInit()); 
    }
    
    delete data;

    return results;  
} 

HpsFitResult* BumpHunter::fitWindow(TH1* histogram, double ap_hypothesis) { 

    // Set the range of the mass variable based on the range of the histogram. 
    variable_map["invariant mass"]->setRange(histogram->GetXaxis()->GetXmin(), histogram->GetXaxis()->GetXmax()); 

    // Create a histogram object compatible with RooFit.
    RooDataHist* data = new RooDataHist("data", "data", RooArgList(*variable_map["invariant mass"]), histogram);

    // Construct a window of size x*(A' mass resolution) around the A' mass
    // hypothesis and do a Poisson likelihood fit within the window range. 
    HpsFitResult* result = this->fitWindow(data, ap_hypothesis);
    
    // Delete the histogram object from memory
    delete data;

    // Return the result
    return result;
}

HpsFitResult* BumpHunter::fitWindow(RooDataHist* data, double ap_hypothesis) { 

    // If the A' hypothesis is below the lower bound, throw an exception.  A 
    // fit cannot be performed using an invalid value for the A' hypothesis.
    this->printDebug("A' mass hypothesis: " + std::to_string(ap_hypothesis)); 
    if (ap_hypothesis < low_bound) throw std::runtime_error("A' hypothesis is less than the lower bound!"); 

    // Get the mass resolution at the mass hypothesis
    double mass_resolution = this->getMassResolution(ap_hypothesis);

    // If the window is being allowed to vary, calculate the window size based
    // on the mass resolution.
    if (!fix_window) { 
        window_size = std::trunc(mass_resolution*15*10000)/10000 + 0.00005;
        
        // If the window size is larger than the max size, set the window size
        // to the max.
        if (window_size > max_window_size) {
            this->printDebug("Window size exceeds maximum."); 
            window_size = max_window_size; 
        }
    }

    // Find the starting position of the window
    double window_start = ap_hypothesis - window_size/2;
    
    // Check that the starting edge of the window is above the boundary.  If not
    // set the starting edge of the window to the boundary.  In this case, the
    // A' hypothesis will not be set to the middle of the window.
    if ((this->low_bound != -9999) && (window_start < this->low_bound)) { 
        this->printDebug("Starting edge of window " + std::to_string(window_start) + " is below lower bound.");
        this->printDebug("Setting edge to lower bound, " + std::to_string(this->low_bound)); 
        window_start = this->low_bound;
    }

    // Check that the end edge of the window is within the high bound.  If not,
    // set the starting edge such that end edge is equal to the high bound.
    if ((this->high_bound != -9999) && ((window_start + window_size) > this->high_bound)) { 
        this->printDebug("End of window " + std::to_string(window_start + window_size) + " is above high bound.");
        this->printDebug("Setting starting edge to " + std::to_string(high_bound - window_size)); 
        window_start = high_bound - window_size;  
    }

    /*
    std::cout << "A' mass: " << ap_hypothesis << std::endl;
    std::cout << "Mass resolution: " << mass_resolution << std::endl;
    std::cout << "Window size: " << window_size << std::endl;
    std::cout << "Window start: " << window_start << std::endl;
    */

    // Set the mean of the Gaussian signal distribution
    variable_map["A' mass"]->setVal(ap_hypothesis);
    
    // Set the width of the Gaussian signal distribution
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

    // Fit the distribution in the given range
    HpsFitResult* result = this->fit(data, false, range_name); 

    // Set the window size 
    result->setWindowSize(window_size);

    // Check if the resulting fit found a significant bump
    double alpha = 0.05; 
    this->calculatePValue(data, result, range_name, alpha); 

    //this->calculateUpperLimit(0.05, data, result, range_name); 

    return result;  
} 

HpsFitResult* BumpHunter::fit(RooDataHist* data, bool migrad_only = false, std::string range_name = "") { 
   
    // Construct a log likelihood using the data set within the range specified
    // above.  This is equivalent to saying that 
    // nll = -ln(L( window_start < x < window_start + window_size | mu, theta))
    // where mu is the signal yield and theta represents the set of all 
    // nuisance parameters which in this case are the background normalization
    // and polynomial constants.  Since the likelihood is being constructed
    // from a histogram, use an extended likelihood.
    //RooAbsReal* nll = model->createNLL(*data, cmd_list);  
    RooAbsReal* nll = model->createNLL(*data, 
            RooFit::Extended(kTRUE), 
            RooFit::Verbose(kTRUE), 
            RooFit::Range(range_name.c_str()), 
            RooFit::SumCoefRange(range_name.c_str())
            );  

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
    /*if (!migrad_only && status == 0) { 
        m.improve();
        m.hesse();
    //    m.minos(*variable_map["signal yield"]); 
    }*/

    // Save the results of the fit
    RooFitResult* result = m.save(); 

    // Delete the constructed negative log likelihood
    delete nll; 

    // Return the saves result
    return new HpsFitResult(result); 
}


void BumpHunter::calculatePValue(RooDataHist* data, HpsFitResult* result, std::string range_name, double alpha) { 

    //  Get the signal yield obtained from the composite fit
    double signal_yield = result->getParameterVal("signal yield");
    
    // We only care if a signal yield is greater than 0.  In the case that it's
    // less than 0, the p-value is set to 1.
    if (signal_yield <= 0) { 
        result->setPValue(1);
        return; 
    }

    // Get the NLL obtained by minimizing the composite model with the signal
    // yield floating.
    double min_nll = result->getRooFitResult()->minNll();

    // Calculate the NLL when signal yield = 0, which is the null hypothesis.
    variable_map["signal yield"]->setVal(0);
    
    // Fix the signal yield at 0.
    variable_map["signal yield"]->setConstant(kTRUE);
   
    // Do the fit
    HpsFitResult* null_result = this->fit(data, true, range_name);

    // Get the NLL obtained from the Bkg only fit.
    double min_nll_null = null_result->getRooFitResult()->minNll(); 
   
    // 1) Calculate the likelihood ratio whose underlying distribution is a 
    //    chi2.
    // 2) From the chi2, calculate the p-value.
    double q0 = 0; 
    double p_value = 0; 
    this->getChi2Prob(min_nll_null, min_nll, q0, p_value);  

    // If P-value is less than the significance, alpha, a bump was found.
    if (p_value < alpha) std::cout << "WTF, a Bump was found!" << std::endl;

    result->setPValue(p_value);
    result->setQ0(q0);  
    
    variable_map["signal yield"]->setConstant(kFALSE);

    delete null_result; 
}

void BumpHunter::printDebug(std::string message) { 
    if (debug) std::cout << "[ BumpHunter ]: " << message << std::endl;
}

void BumpHunter::resetParameters(RooArgList initial_params) { 
    
    for (auto& element : variable_map) { 
        if (initial_params.find(element.second->GetName()) == NULL) continue;

        RooRealVar* var = (RooRealVar*) initial_params.find(element.second->GetName());

        element.second->setVal(var->getVal());
        element.second->setError(var->getError()); 

    }
}

/*void BumpHunter::calculateUpperLimit(double alpha, RooDataHist* data, HpsFitResult* result, std::string range_name) { 

    
    // Get the value of the NLL at the minimum.
    // Now search for the upper limit
    while(true) {

        // Set the signal yield constant.
        variable_map["signal yield"]->setConstant(kFALSE);
    
        // Calculate the NLL for the signal yield = 0, which is the null hypothesis.
        std::cout << "Prob: " << prob << std::endl;
        if (prob > 0.5) signal_yield += 50;
        else if (prob > 0.2) signal_yield += 10;
        else if (prob > 0.1) signal_yield++;   
        //std::cout << "Signal yield: " << signal_yield << std::endl;
        //signal_yield++; 
        //variable_map["signal yield"]->setVal(signal_yield);
        //std::cout << "New signal yield: " << signal_yield << std::endl;

        // Set the signal yield constant.
       variable_map["signal yield"]->setConstant(kTRUE);

        HpsFitResult* current_result = this->fit(data, true, range_name);
        
        min_nll = current_result->getRooFitResult()->minNll(); 

        prob = this->getChi2Prob(min_nll_null, min_nll);  
        if (prob < 0.1) { 
            std::cout << "Upper limit: " << signal_yield << std::endl;
            result->setUpperLimit(signal_yield); 
            delete current_result; 
            break; 
        }
        delete current_result; 
    }

    //delete null_result;

    // Set the signal yield constant.
    variable_map["signal yield"]->setConstant(kFALSE);
    
}*/

void BumpHunter::getChi2Prob(double min_nll_null, double min_nll, double &q0, double &p_value) {
    //std::cout << std::fixed << "Null NLL: " << min_nll_null << std::endl;
    //std::cout << std::fixed << "Min NLL: " << min_nll << std::endl;
    
    min_nll *= -1;
    min_nll_null *= -1;

    double diff = min_nll - min_nll_null;
    //std::cout << "Difference: " << diff << std::endl;
    
    q0 = 2*diff;
    //std::cout << "q0: " << q0 << std::endl;
    
    p_value = TMath::Prob(q0, 1); 
    //std::cout << "Probability: " << prob << std::endl;
}

void BumpHunter::fitBkgOnly() { 
    this->bkg_only = true; 
    model = bkg_model; 
}

void BumpHunter::setBounds(double low_bound, double high_bound) {
    this->low_bound = low_bound; 
    this->high_bound = high_bound;
    printf("Fit bounds set to [ %f , %f ]\n", this->low_bound, this->high_bound);   
}

void BumpHunter::writeResults() { 
    
    // Create the output file name string
    char buffer[100];
    sprintf(buffer, "results_order%i_window%i.txt", bkg_poly_order, window_size*1000);

    // Create a file stream  
    ofs = new std::ofstream(buffer, std::ofstream::out); 
}
