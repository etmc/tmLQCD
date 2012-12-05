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

/******************************************
 * Hopping_Matrix is the conventional Wilson 
 * hopping matrix
 *
 *
 * But the communication is left out by a 
 * dirty trick...
 *
 * \kappa\sum_{\pm\mu}(r+\gamma_\mu)U_{x,\mu}
 *
 * for ieo = 0 this is M_{eo}, for ieo = 1
 * it is M_{oe}
 *
 * l is the number of the output field
 * k is the number of the input field
 *
 ******************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "xchange/xchange.h"
#include "su3.h"
#include "sse.h"
#include "boundary.h"
#include "operator/Hopping_Matrix.h"

#define Hopping_Matrix Hopping_Matrix_nocom
#define _NO_COMM 1
#ifdef _KOJAK_INST
#undef _KOJAK_INST
#endif

#include "Hopping_Matrix.c"
