/**
 *
 *
 *
 */

#ifndef __HPS_FIT_RESULT_H__
#define __HPS_FIT_RESULT_H__

//------------//
//   RooFit   //
//------------//
#include <RooFitResult.h>
#include <RooRealVar.h>

class HpsFitResult { 

    public: 

        /** Default constructor */
        HpsFitResult(); 

        /** */
        HpsFitResult(RooFitResult* result, double q0 = 0, double p_value = 0, double upper_limit = 0);

        ~HpsFitResult(); 

        /** */
        double getQ0() { return q0; };

        /** */
        double getPValue() { return p_value; };

        /** */
        double getParameterVal(std::string parameter_name); 

        /** */
        RooFitResult* getRooFitResult() { return result; };  

        double getBkgTotal() {return this->bkg_total; };

        double getBkgWindowSize() { return this->bkg_window_size; };  


        /** */
        double getUpperLimit() { return upper_limit; };

        /** @return The size of the fit window used. */
        double getWindowSize() { return this->window_size; }; 
        
        /** */
        double setQ0(double q0) { this->q0 = q0; };

        /** */
        void setPValue(double p_value) { this->p_value = p_value; };  
        
        /** */
        void setRooFitResult(RooFitResult* result) { this->result = result; }; 

        void setBkgTotal(double bkg_total) { this->bkg_total = bkg_total; };

        void setBkgWindowSize(double bkg_window_size) { this->bkg_window_size = bkg_window_size; };  

        /**
         * Set the 2 sigma upper limit.
         *
         * @param upper_limit The 2 sigma upper limit.
         */
        void setUpperLimit(double upper_limit) { this->upper_limit = upper_limit; };

        /**
         * Set the size of the fit window used to get this results.
         *
         * @oaram window_size The size of the fit window.
         */
        void setWindowSize(double window_size) { this->window_size = window_size; }; 

    private: 

        /** Result associated with RooFit. */
        RooFitResult* result; 

        double bkg_total;

        double bkg_window_size; 

        /** q0 value */
        double q0;
        
        /** p-value. */
        double p_value;

        /** 2 sigma upper limit on the signal. */
        double upper_limit; 

        /*** Size of the fit window. */
        double window_size;

}; // HpsFitResult

#endif // __HPS_FIT_RESULT_H__
