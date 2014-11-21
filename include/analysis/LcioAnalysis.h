/**
 *	@file LcioAnalysis.h
 *	@brief Interface describing an LCIO analysis.
 *	@author Omar Moreno <omoreno1@ucsc.edu>
 *	@date November 17, 2014
 *
 */

#ifndef __LCIO_ANALYSIS_H__
#define __LCIO_ANALYSIS_H__

#include <string>

#include <EVENT/LCEvent.h>

class LcioAnalysis {

	public: 

		virtual ~LcioAnalysis() { }; 

		/**
		 *
		 */
		virtual void initialize() = 0; 

		/**
		 *
		 */
		virtual void processEvent(EVENT::LCEvent*) = 0;

		/**
		 *
		 */
		virtual void finalize() = 0; 

		/**
		 *
		 */
		virtual std::string toString() = 0; 

}; // LcioAnalysis 

#endif // __LCIO_ANALYSIS_H__
