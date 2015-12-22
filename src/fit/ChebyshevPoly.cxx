/**
 *
 */

#include <ChebyshevPoly.h>

double ChebyshevPoly::firstOrder(double x, double par[]) {
    return par[0] + par[1] * x; 
}



double ChebyshevPoly::firstOrderIntegral(double x, double par[]) {   
   return par[0]*x + (par[1]*x*x)/2; 
}

double ChebyshevPoly::firstOrderIntegral(double x_low, double x_high, double par[]) { 
    return firstOrderIntegral(x_high, par) - firstOrderIntegral(x_low, par); 
}
