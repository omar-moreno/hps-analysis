/**
 *
 */

#include <ChebyshevPoly.h>

double ChebyshevPoly::firstOrder(double x, double par[]) {
    return par[0] + par[1] * x; 
}
