
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
    RooDataHist* data = new RooDataHist("data", "data", RooArgList(*variable_map["invariant mass"]), histogram);

    // Create a container for the results from the fit to each window.
    std::map<double, RooFitResult*> results; 

    while (start <= (end - window_size)) { 
        
        double ap_hypothesis = start + window_size/2; 
        variable_map["A' mass"]->setVal(ap_hypothesis);
        variable_map["A' mass resolution"]->setVal(this->getMassResolution(ap_hypothesis)); 

        std::string range_name = "ap_mass_" + std::to_string(ap_hypothesis); 
        variable_map["invariant mass"]->setRange(range_name.c_str(), start, start + window_size); 

        double min_bin = histogram->GetXaxis()->FindBin(start); 
        double max_bin = histogram->GetXaxis()->FindBin(start + window_size);
        double integral = histogram->Integral(min_bin, max_bin);  
        std::cout << "Estimated bkg in range (" << start << ", " << start + window_size << "): " << integral << std::endl;
        variable_map["bkg yield"]->setVal(integral);

        RooAbsReal* nll = model->createNLL(*data, 
                RooFit::Extended(kTRUE), 
                RooFit::SumCoefRange(range_name.c_str()), 
                RooFit::Range(range_name.c_str()),
                RooFit::Verbose(kFALSE));

        RooMinuit m(*nll);

        m.setPrintLevel(-1000);

        // Set the Minuit strategy
        m.setStrategy(1);

        if (m.migrad() != 0) { 
            m.simplex();
            m.migrad();
        }

        m.improve();

        m.hesse();

        m.minos(*variable_map["signal yield"]);

        RooFitResult* result = m.save(); 
        results[ap_hypothesis] = result; 
        if (ofs != nullptr) result->printMultiline(*ofs, 0, kTRUE, "");

        start += window_step; 

        resetParameters(result->floatParsInit()); 
         
        delete nll; 
    }

    delete data; 
    return results; 
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
