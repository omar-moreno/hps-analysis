/**
 *	@author:	Omar Moreno <omoreno1@ucsc.edu>
 *	@section institution
 *				Santa Cruz Institute for Particle Physics
 *				University of California, Santa Cruz
 *	@date:		March 31, 2014 
 *
 */

#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__

class Topology { 

	public: 

		virtual ~Topology(){}; 

		virtual void setup() = 0; 
		virtual void cuts() = 0; 
		virtual void analyze(HpsEvent*) = 0; 
		virtual void finalize() = 0;

}; // Topology

#endif // __TOPOLOGY_H__
