/***************************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2017                               Bartosz Kostrzewa
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
 ****************************************************************************/

#ifndef _SOLVER_TYPES_H
#define _SOLVER_TYPES_H

typedef enum SOLVER_TYPE {
 BICGSTAB = 0,
 CG,
 GMRES,
 CGS,
 MR,
 BICGSTABELL,
 FGMRES,
 GCR,
 GMRESDR,
 PCG,
 DFLGCR,
 DFLFGMRES,
 CGMMS,
 MIXEDCG,
 RGMIXEDCG,
 CGMMSND,
 INCREIGCG,
 MIXEDCGMMSND,
 SUMR,
 MCR,
 CR,
 BICG,
 MG,
 MIXEDBICGSTAB,
 DUMMYHERMTEST
} SOLVER_TYPE;

int solver_is_mixed( const int solver_type );

#endif
