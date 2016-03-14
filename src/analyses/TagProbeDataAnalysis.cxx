/**
 * @file TagProbeDataAnalysis.cxx
 * @brief Analysis used to calculate tracking efficiency using Moller events
 *        in the Engineering Run 2015 data.
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date February 16, 2015
 */

#include <TagProbeDataAnalysis.h>

void TagProbeDataAnalysis::processEvent(HpsEvent* event) { 
    
    // Only look at pairs1 triggers
    if (!event->isSingle1Trigger()) return;

    // Only look at events with the SVT bias ON
    if (!event->isSvtBiasOn()) return; 

    // Only look at events where the SVT is closed
    if (!event->isSvtClosed()) return;

    // Skip events that had burst mode noise
    if (event->hasSvtBurstModeNoise()) return; 

    // Skip events that had SVT header errors
    if (event->hasSvtEventHeaderErrors()) return;

    // 
    TagProbeAnalysis::processEvent(event); 
}
