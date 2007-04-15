/* $Id$ */

/****************************************************************
 *
 * invert_eo makes an inversion with EO precoditioned
 * tm Operator
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ****************************************************************/

#ifndef _INVERT_EO_H
#define _INVERT_EO_H

int invert_eo(spinor * const Even_new, spinor * const Odd_new, 
	      spinor * const Even, spinor * const Odd,
	      const double precision, const int iter_max,
	      const int solver_flag, const int rel_prec,
	      const int sub_evs_flag);

/* This is the full matrix multiplication */
void M_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd);
void Q_full(spinor * const Even_new, spinor * const Odd_new, 
	    spinor * const Even, spinor * const Odd);
void M_minus_1_timesC(spinor * const Even_new, spinor * const Odd_new, 
		      spinor * const Even, spinor * const Odd);
#endif
