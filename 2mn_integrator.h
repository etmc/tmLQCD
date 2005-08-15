/* $Id$ */
#ifndef _2MN_INTEGRATOR_H
#define _2MN_INTEGRATOR_H

/**************************************************
 * 
 * velocity version 
 *
 **************************************************/

void mn2_integrator(int * const n_int, const double tau, 
		    const int S, const int halfstep, double * const lambda);

/**************************************************
 * 
 * position version 
 *
 **************************************************/


void mn2p_integrator(int * const n_int, const double tau, 
		     const int S, double * const lambda);
#endif
