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

std::vector<double> f_vals; 

std::vector<double> t0;

std::vector<double> t1; 

/** Bin width */
double bin_width = 0; 

double histo_integral = 0; 

void getData(TH1* histo) { 
    bin_width = histo->GetXaxis()->GetBinWidth(1);
    histo_integral = histo->GetEntries(); 
    //std::cout << "Histogram integral: " << histo_integral << std::endl;
    for (int bin_n = 1; bin_n < histo->GetXaxis()->GetNbins(); ++bin_n) { 
        x_vec[histo->GetXaxis()->GetBinCenter(bin_n)] = histo->GetBinContent(bin_n); 
    }
}


void fcn(int& npar, double* deriv, double& f, double par[], int flag) {
    double lnL = 0.0; 

    double norm = ChebyshevPoly::firstOrderIntegral(0, 10, par);
    //std::cout << "Normalization: " << norm << std::endl; 
    for (auto& x : x_vec) {
        double x_low  = x.first - bin_width/2;
        double x_high = x.first + bin_width/2;
        //std::cout << "x low: "  << x_low << std::endl;
        //std::cout << "x high: " << x_high << std::endl;
        double func_value = histo_integral*(ChebyshevPoly::firstOrderIntegral(x_low, x_high, par)/norm); 
        //std::cout << "v(theta) = " << func_value << std::endl; 

        if (func_value > 0) {
            lnL += x.second*log(func_value);
        } else { 
            //std::cout << "WARNING! Value is negative." << std::endl;
        }
    }

    f = -lnL + histo_integral; 
    
    f_vals.push_back(f);
    t0.push_back(par[0]); 
    t1.push_back(par[1]); 
    
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
    chebyshev->SetParameter(0, 10);
    chebyshev->SetParameter(1, -1);

    std::cout << "Integral: " << chebyshev->Integral(0, 10) << std::endl;

    TH1F* toy_histo = new TH1F("toy_histo", "toy_histo", 100, 0, 10); 
    toy_histo->FillRandom("chebyshev", 100000);
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
    
    par[0] = 15.0; 
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

    minuit->mnhess();

    //minuit->mnimpr(); 

    TF1* fit_result = new TF1("fit_result", "(100000/50*((10 - x)/50))", 0, 10); 

    toy_histo->Draw();
    chebyshev->Draw("same");
    fit_result->Draw("same"); 
    
    canvas->SaveAs("chebyshev_fit_example.png");

    TH1F* f_histo = new TH1F("lnL", "lnL", 50, -611000, -602000); 
    for (auto f_val : f_vals) { 
        f_histo->Fill(f_val); 
    }

    f_histo->Draw(); 

    canvas->SaveAs("lnL.png");

    TH1F* t0_histo = new TH1F("t0", "t0", 50, 9.8, 10.4); 
    for (auto t0_val : t0) { 
        t0_histo->Fill(t0_val); 
    }

    t0_histo->Draw(); 

    canvas->SaveAs("t0.png");

    TH1F* t1_histo = new TH1F("t1", "t1", 50, -20, 20); 
    for (auto t1_val : t1) { 
        t1_histo->Fill(t1_val); 
    }

    t1_histo->Draw(); 

    canvas->SaveAs("t1.png");



    delete toy_histo;
    delete chebyshev; 
    delete canvas; 
}


