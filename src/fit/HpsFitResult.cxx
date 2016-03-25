
#include <HpsFitResult.h>

HpsFitResult::HpsFitResult() 
    : result(nullptr), 
      upper_limit(0), 
      p_value(0) { 
}

HpsFitResult::HpsFitResult(RooFitResult* result, double upper_limit, double p_value)  
    : result(result), 
      upper_limit(upper_limit), 
      p_value(p_value) { 
}

HpsFitResult::~HpsFitResult() { 
    if (result != nullptr) delete result;
}

double HpsFitResult::getParameterVal(std::string parameter_name) { 
    return ((RooRealVar*) this->getRooFitResult()->floatParsFinal().find(parameter_name.c_str()))->getVal(); 
}
