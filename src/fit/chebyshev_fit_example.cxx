/**
 *
 */

#include <map>

//----------//
//   ROOT   //
//----------//
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMinuit.h>

#include <ChebyshevPoly.h>

/** */
std::map<double, double> x_vec;

/** Bin width */
double bin_width = 0; 

double histo_integral = 0; 

void getData(TH1* histo) { 
    bin_width = histo->GetXaxis()->GetBinWidth(1);
    histo_integral = histo->GetEntries(); 
    for (int bin_n = 1; bin_n < histo->GetXaxis()->GetNbins(); ++bin_n) { 
        x_vec[histo->GetXaxis()->GetBinCenter(bin_n)] = histo->GetBinContent(bin_n); 
    }
}


void fcn(int& npar, double* deriv, double& f, double par[], int flag) {
    double lnL = 0.0; 

    double total = ChebyshevPoly::firstOrder(0, par) - ChebyshevPoly::firstOrder(10, par); 
    for (auto& x : x_vec) {
        double x_low  = x.first - bin_width/2;
        double x_high = x.first + bin_width/2;
        double func_value =  histo_integral*(ChebyshevPoly::firstOrder(x_low, par) 
                                             - ChebyshevPoly::firstOrder(x_high, par));

        if (func_value > 0) {
            lnL += x.second*log(func_value) - total;
            //lnL += log(func_value);
        } else { 
            std::cout << "WARNING! Value is negative." << std::endl;
        }
    }

    f = -lnL;
    //std::cout << "f: " << f << std::endl; 
}

void setupCanvas(TCanvas* canvas) {
    canvas->SetFillColor(0);
 	canvas->SetBorderMode(0);
 	canvas->SetBorderSize(0);
 	canvas->SetFrameFillColor(0);
    canvas->SetFrameBorderMode(0);
}

int main(int argc, char** argv) {
  
    TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
    setupCanvas(canvas);
    
    TF1* chebyshev = new TF1("chebyshev", "[0] + [1]*x", 0, 10);
    chebyshev->SetParameter(0, 20);
    chebyshev->SetParameter(1, -1);

    TH1F* toy_histo = new TH1F("toy_histo", "toy_histo", 100, 0, 10); 
    toy_histo->FillRandom("chebyshev", 1000);
    toy_histo->GetXaxis()->SetTitle("x");
    toy_histo->GetYaxis()->SetTitle("Events/100");
   
    getData(toy_histo); 

    //
    // Minuit
    //
    
    const int npar = 2; 
    
    TMinuit* minuit = new TMinuit(npar);
    minuit->SetPrintLevel(1);  
    minuit->SetFCN(fcn);

    std::cout << "Iterations: " << minuit->GetMaxIterations() << std::endl;

    double par[npar]; 
    double step_size[npar]; 
    double min_val[npar]; 
    double max_val[npar];
    std::string par_name[npar]; 
    
    par[0] = 40.0; 
    step_size[0] = 1; 
    min_val[0] = 0;
    max_val[0] = 0; 
    par_name[0] = "t0"; 

    par[1] = -1; 
    step_size[1] = .2; 
    min_val[1] = 0;
    max_val[1] = 0; 
    par_name[1] = "t1";

    for (int par_n = 0; par_n < npar; ++par_n) { 
        minuit->DefineParameter(par_n, par_name[par_n].c_str(), 
                par[par_n], step_size[par_n], min_val[par_n], 
                max_val[par_n]); 
    } 

    minuit->Migrad();

    toy_histo->Draw();
    chebyshev->Draw("same");
    
    canvas->SaveAs("chebyshev_fit_example.png"); 
    
    delete toy_histo;
    delete chebyshev; 
    delete canvas; 
}


