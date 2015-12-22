/**
 *
 */

#ifndef __CHEBYSHEV_POLY_H__
#define __CHEBYSHEV_POLY_H__

namespace ChebyshevPoly { 

    /**
     *
     */
    double firstOrder(double x, double par[]);

    /**
     *
     */
    double firstOrderIntegral(double x, double par[]);

    /**
     *
     */
    double firstOrderIntegral(double x_low, double x_high, double par[]); 
}

#endif // __CHEBYSHEV_POLY_H__
