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
/* 
 *  Routine for the computation of the Jacobi operator (for use into LapH_ev)
 *  Authors Luigi Scorzato, Marco Cristoforetti
 *
 *
 *******************************************************************************/
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "su3.h"
#include "xchange/xchange.h"

#ifdef WITHLAPH

void Jacobi(su3_vector * const l, su3_vector * const k,int t)
{
  int ix,mu,tcoord,coord;
  su3_vector lt;
        
#ifdef MPI
  xchange_jacobi(k);
#endif

  tcoord=t*SPACEVOLUME;
  for(ix=0;ix<SPACEVOLUME;ix++)
    {
      coord=tcoord+ix;
      _vector_mul(l[ix],6,k[ix]);
      for(mu=1;mu<4;mu++)
	{
	  _su3_multiply(lt,g_gauge_field[coord][mu],k[g_iup3d[ix][mu]]);
	  l[ix].c0 -= lt.c0;
	  l[ix].c1 -= lt.c1;
	  l[ix].c2 -= lt.c2;
	  _su3_inverse_multiply(lt,g_gauge_field[g_idn[coord][mu]][mu],k[g_idn3d[ix][mu]]);
	  l[ix].c0 -= lt.c0;
	  l[ix].c1 -= lt.c1;
	  l[ix].c2 -= lt.c2;
	}
    }
#ifdef MPI
  xchange_jacobi(l);
#endif
}

#endif // WITHLAPAH
