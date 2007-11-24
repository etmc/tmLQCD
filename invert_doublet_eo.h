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

#ifndef _INVERT_DOUBLET_EO_H
#define _INVERT_DOUBLET_EO_H

int invert_doublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
		      spinor * const Even_new_c, spinor * const Odd_new_c, 
		      spinor * const Even_s, spinor * const Odd_s,
		      spinor * const Even_c, spinor * const Odd_c,
		      const double precision, const int max_iter,
		      const int solver_flag, const int rel_prec);


/* This is the full matrix multiplication */
/* void M_full(spinor * const Even_new, spinor * const Odd_new,  */
/* 	    spinor * const Even, spinor * const Odd); */
/* void Q_full(spinor * const Even_new, spinor * const Odd_new,  */
/* 	    spinor * const Even, spinor * const Odd); */
/* void M_minus_1_timesC(spinor * const Even_new, spinor * const Odd_new,  */
/* 		      spinor * const Even, spinor * const Odd); */
#endif
