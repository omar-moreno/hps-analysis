
#ifndef __V0_DATA_ANALYSIS_H__
#define __V0_DATA_ANALYSIS_H__

#include <V0Analysis.h>

class V0DataAnalysis : public V0Analysis { 

    void processEvent(HpsEvent* event); 

    std::string toString(); 

};

#endif
