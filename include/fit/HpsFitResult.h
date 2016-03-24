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
        HpsFitResult(RooFitResult* result, double upper_limit = 0, double p_value = 0); 

        ~HpsFitResult(); 

        double getPValue() { return p_value; };

        /** 
         *
         */
        double getParameterVal(std::string parameter_name); 


        RooFitResult* getRooFitResult() { return result; };  

        double getUpperLimit() { return upper_limit; };
        
        /**
         * Set the 2 sigma upper limit.
         *
         * @param upper_limit The 2 sigma upper limit.
         */
        void setUpperLimit(double upper_limit) { this->upper_limit = upper_limit; };

        /**
         *
         */
        void setPValue(double p_value) { this->p_value = p_value; };  

        void setRooFitResult(RooFitResult* result) { this->result = result; }; 

    private: 

        /** Result associated with RooFit. */
        RooFitResult* result; 

        /** 2 sigma upper limit on the signal. */
        double upper_limit; 

        /** p value. */
        double p_value;

}; // HpsFitResult

#endif // __HPS_FIT_RESULT_H__
