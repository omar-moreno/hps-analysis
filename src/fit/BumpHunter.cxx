
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
    : window_size(0.01) {

    // Independent variable
    variable_map["invariant mass"] = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0.003, 0.1);

    //   Signal PDF   //
    variable_map["A' mass"]  = new RooRealVar("A' Mass",  "A' Mass",  0.03);
    variable_map["A' mass resolution"] = new RooRealVar("A' Mass Resolution", "A' Mass Resolution", 0.00167);

    signal = new RooGaussian("signal", "signal", *variable_map["invariant mass"],
                             *variable_map["A' mass"], *variable_map["A' mass resolution"]);

    std::string name;
    for (int order = 1; order <= poly_order; ++order) {
        name = "t" + std::to_string(order);
        variable_map[name] = new RooRealVar(name.c_str(), name.c_str(), 0, -1, 1);
        arg_list.add(*variable_map[name]);
    }

    
    bkg = new RooChebychev("bkg", "bkg", *variable_map["invariant mass"], arg_list);


    variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -5000, 5000);
    variable_map["bkg yield"] = new RooRealVar("bkg yield", "bkg yield", 50000, 10000, 1000000);

    model = new RooAddPdf("model", "model", RooArgList(*signal, *bkg), 
                            RooArgList(*variable_map["signal yield"], *variable_map["bkg yield"])); 
}


BumpHunter::~BumpHunter() {

    
    for (auto& element : variable_map) { 
       delete element.second; 
    }
    variable_map.clear();

    delete signal;
    delete bkg;
    delete model; 
}
