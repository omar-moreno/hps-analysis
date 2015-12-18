/**
 *
 */

//----------//
//   ROOT   //
//----------//
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>

/*vector<double> x_vec; 

void fcn(int& npar, double* deriv, double& f, double par[], int flag) {

}*/

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
    
    TF1* chebyshev = new TF1("chebyshev", "[0] + [1]*x", -10, 10);
    chebyshev->SetParameter(0, 1);
    chebyshev->SetParameter(1, 5);

    TH1F* toy_histo = new TH1F("toy_histo", "toy_histo", 100, 0, 10); 
    for (int event_n = 0; event_n < 10000; ++event_n) { 
        toy_histo->Fill(chebyshev->GetRandom(0, 10));
    }
     
    toy_histo->Draw();
    //chebyshev->Draw("same");
    
    canvas->SaveAs("chebyshev_fit_example.png"); 
    
    delete toy_histo;
    delete chebyshev; 
    delete canvas; 
}


