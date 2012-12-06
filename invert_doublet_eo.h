/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

/****************************************************************
 *
 * invert_doublet_eo makes an inversion with EO precoditioned
 * tm Operator with a nondegenerate doublet
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

int invert_cloverdoublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
			    spinor * const Even_new_c, spinor * const Odd_new_c, 
			    spinor * const Even_s, spinor * const Odd_s,
			    spinor * const Even_c, spinor * const Odd_c,
			    const double precision, const int max_iter,
			    const int solver_flag, const int rel_prec);
#endif
