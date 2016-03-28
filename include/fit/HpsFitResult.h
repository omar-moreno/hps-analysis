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

        /** */
        double getUpperLimit() { return upper_limit; };
        
        /** */
        double setQ0(double q0) { this->q0 = q0; };

        /** */
        void setPValue(double p_value) { this->p_value = p_value; };  
        
        /** */
        void setRooFitResult(RooFitResult* result) { this->result = result; }; 

        /**
         * Set the 2 sigma upper limit.
         *
         * @param upper_limit The 2 sigma upper limit.
         */
        void setUpperLimit(double upper_limit) { this->upper_limit = upper_limit; };


    private: 

        /** Result associated with RooFit. */
        RooFitResult* result; 

        /** q0 value */
        double q0;
        
        /** p value. */
        double p_value;

        /** 2 sigma upper limit on the signal. */
        double upper_limit; 

}; // HpsFitResult

#endif // __HPS_FIT_RESULT_H__