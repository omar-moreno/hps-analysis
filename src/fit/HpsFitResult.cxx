
#include <HpsFitResult.h>

HpsFitResult::HpsFitResult() 
    : result(nullptr),
      q0(0),  
      p_value(0),
      upper_limit(0) { 
}

HpsFitResult::HpsFitResult(RooFitResult* result, double q0, double p_value, double upper_limit)  
    : result(result), 
      q0(q0),
      p_value(p_value),
      upper_limit(upper_limit) { 
}

HpsFitResult::~HpsFitResult() { 
    if (result != nullptr) delete result;
}

double HpsFitResult::getParameterVal(std::string parameter_name) { 
    return ((RooRealVar*) this->getRooFitResult()->floatParsFinal().find(parameter_name.c_str()))->getVal(); 
}
