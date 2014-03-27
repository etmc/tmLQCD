/***********************************************************************
 * Copyright (C) 2014 Florian Burger
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
 *
 * invert_noeo makes an inversion with non-EO precoditioned
 * tm Operator
 *
 * Author: Florian Burger
 *         burger@physik.hu-berlin.de
 *
 ***********************************************************************/

#ifndef _INVERT_NOEO_H
#define _INVERT_NOEO_H

int invert_noeo(spinor * const Spin_new,  
	      spinor * const Spin, 
	      const double precision, const int iter_max,
	      const int solver_flag, const int rel_prec,
              const int id );

#endif
