'''

'''

###############
#   Imports   #
###############
import ROOT as r

class BumpHunter :

    def __init__(self, poly_order) : 

        # Create the object used to define the observable
        self.invariant_mass = r.RooRealVar("Invariant Mass",
                                           "Invariant Mass (GeV)", 0.005, 0.1) 
        # Define the window size that will be used in the fit
        self.mass_window_size = 0.01

        #
        # Build the composite model
        #
        
        #
        # Signal PDF
        #

        # Create the objects that represent the mass hypothesis and mass
        # resolution. These will be used in defining the Gaussian model used to
        # descibe the A' signal.  These values are set constant during the fit.
        self.ap_mass_mean = r.RooRealVar("ap_mass_mean", "ap_mass_mean", .03)
        #self.ap_mass_mean.setConstant(r.kTRUE)
        self.ap_mass_sigma = r.RooRealVar("ap_mass_sigma", "ap_mass_sigma", 0.00167)
        #self.ap_mass_sigma.setConstant(r.kTRUE)

        # Define the Gaussian model used to describe the A' signal.
        self.signal = r.RooGaussian("signal", "signal",
                                    self.invariant_mass, 
                                    self.ap_mass_mean,
                                    self.ap_mass_sigma)

        # Create the objects that represent the polynomial constants used to
        # define the background. These variables are stored in a list for easy
        # access.
        self.t = []
        self.arg_list = r.RooArgList()

        for order in range(1, poly_order+1) :
            self.t.append(r.RooRealVar("t"+str(order), "t"+str(order), 0, -10000, 10000))
            self.arg_list.add(self.t[order - 1])

        # Define the polynomial model used to describe the background in some
        # pre-determined window.
        self.bkg = r.RooChebychev("bkg", "bkg", self.invariant_mass, self.arg_list)

        # Create the objects that will represent the number of signal and 
        # background events in a given window.
        self.nsig = r.RooRealVar("nsig","signal fraction", 0, -1000., 1000)
        self.nbkg = r.RooRealVar("nbkg","background fraction", 10000., 0.,1000000000)

        # Build a composite model 
        self.model = r.RooAddPdf("model", "model", 
                r.RooArgList(self.signal, self.bkg), r.RooArgList(self.nsig, self.nbkg))

    
    def fit(self, histogram, mass_start, mass_end, mass_step) :

        histogram_data = r.RooDataHist("invariant_mass_data", "invariant_mass_data", 
                r.RooArgList(self.invariant_mass), histogram)

        # List used to store all of the fit results
        results = []

        # Loop over 
        while mass_start <= mass_end - self.mass_window_size : 

            ap_mass = mass_start + self.mass_window_size/2
            self.ap_mass_mean.setVal(ap_mass)
            
            self.invariant_mass.setRange("A' mass = " + str(ap_mass),
                                         mass_start, mass_start + self.mass_window_size)

            nll = self.model.createNLL(histogram_data, 
                                       r.RooFit.Extended(r.kTRUE), 
                                       r.RooFit.SumCoefRange("A' mass = " + str(ap_mass)),
                                       r.RooFit.Range("A' mass = " + str(ap_mass)))
            m = r.RooMinuit(nll)

            m.migrad()

            m.improve()

            m.hesse()

            #m.minos()

            result = m.save()
            #result = self.model.fitTo(histogram_data, 
            #                           r.RooFit.SumCoefRange("A' mass = " + str(ap_mass)),
            #                           r.RooFit.Range("A' mass = " + str(ap_mass)),
            #                          r.RooFit.Extended(r.kTRUE),
            #                           r.RooFit.Save())

            results.append(result)
            self.reset_params(result.floatParsInit())

            mass_start += mass_step
    
        return results

    
    def reset_params(self, initial_params) : 

        for t_var in self.t : 
        
            #print "Resetting value of " + str(a_var.GetName())
            #print "Current value: " + str(a_var.getVal()) + " Current Error: " + str(a_var.getError())
        
            if initial_params.index(t_var.GetName()) == -1 : 
                #print "Variable is not in argument list."
                continue 
            t_var.setVal(initial_params[initial_params.index(t_var.GetName())].getVal())
            t_var.setError(initial_params[initial_params.index(t_var.GetName())].getError())
        
        self.nsig.setVal(initial_params[initial_params.index(self.nsig.GetName())].getVal())
        self.nsig.setError(initial_params[initial_params.index(self.nsig.GetName())].getError())
        self.nbkg.setVal(initial_params[initial_params.index(self.nbkg.GetName())].getVal())
        self.nbkg.setError(initial_params[initial_params.index(self.nbkg.GetName())].getError())
            #print "Reset value: " + str(a_var.getVal()) + " Reset Error: " + str(a_var.getError())
        

    def set_mass_window_size(self, mass_window_size) : 
        self.mass_window_size = mass_window_size
