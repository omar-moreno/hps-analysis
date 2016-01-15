
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
    : window_size(0.01)
{

    // 
    invariant_mass = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0.003, 0.1);

    //   Signal PDF   //
    ap_mass_mean  = new RooRealVar("ap_mass_mean",  "ap_mass_mean",  0.03);
    ap_mass_sigma = new RooRealVar("ap_mass_sigma", "ap_mass_sigma", 0.00167);

    signal = new RooGaussian("signal", "signal", *invariant_mass, *ap_mass_mean, *ap_mass_sigma);

    for (int order = 1; order <= poly_order; ++order) {
        t.push_back(new RooRealVar(("t"+std::to_string(order)).c_str(),
                    ("t"+std::to_string(order)).c_str(), 
                    0, -1, 1));
        arg_list.add(*t[order - 1]);
    }

    bkg = new RooChebychev("bkg", "bkg", *invariant_mass, arg_list);

    n_sig = new RooRealVar("nsig", "signal yield", 0, -5000, 5000);
    n_bkg = new RooRealVar("nbkg", "bkg yield", 50000, 10000, 1000000);

    model = new RooAddPdf("model", "model", RooArgList(*signal, *bkg), RooArgList(*n_sig, *n_bkg)); 
}


BumpHunter::~BumpHunter() {
    delete invariant_mass;
    delete ap_mass_mean;
    delete ap_mass_sigma;
    delete signal;
}
