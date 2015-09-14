/**
 *
 * @file StrategyType.h
 * @brief Enum constants denoting the different track finding algorithms.
 * @author <a href="mailto:omoreno1@ucsc.edu">Omar Moreno</a>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date September 13, 2015
 *
 */

#ifndef __STRATEGY_TYPE_H__
#define __STRATEGY_TYPE_H__

enum class StrategyType { 
    
    /** Enum values denoting the different track finding algorithsm. */
    /** 
     * This is the default StrategyType value set by the reconstruction.
     * For pass 2, it represents the following strategy:
     *
     * Seed 345, confirm 2, extend 16
     *
     * For pass 3 and beyond, this should no longer be used.
     */
    MATCHED_TRACKS = 1, 
    /** Seed 345, confirm 2, extend 16 */
    S345_C2_E16 = 2, 
    /** Seed 456, confirm 3, extend 21 */
    S456_C3_E21 = 3,
    /** Seed 123, confirm 4, extend 56 */
    S123_C4_E56 = 4, 
    /** Seed 123, confirm 5, extend 46 */
    S123_C5_E45 = 5

}; // Strategy Type

#endif // __STRATEGY_TYPE_H__
