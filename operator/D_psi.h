/***********************************************************************
 *
 * Copyright (C) 2007,2008 Carsten Urbach
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

#ifndef _D_PSI_H
#define _D_PSI_H

#include "block.h"

void D_psi(spinor * const P, spinor * const Q);
void D_psi_prec(spinor * const P, spinor * const Q);
void Block_D_psi(block * blk, spinor * const rr, spinor * const s);
void Block_H_psi(block * blk, spinor * const rr, spinor * const s, const int eo);


#if (defined BGL && defined XLC)
void local_H(spinor * const rr, spinor * const s, su3 * u, int * _idx);
#else
inline void local_H(spinor * const rr, spinor const * const s, su3 const * restrict u, int * _idx, spinor * const restrict tmpr);
#endif

inline void p0add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void m0add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);


inline void p1add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void m1add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void p2add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void m2add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void p3add(spinor * restrict const tmpr , spinor const * restrict const s,
                         su3 const * restrict const u, const _Complex double phase);

inline void m3addandstore(spinor * restrict const r, spinor const * restrict const s,       
                                 su3 const * restrict const u, const _Complex double phase,
         spinor const * restrict const tmpr);



void boundary_D_0(spinor * const r, spinor * const s, su3 *u);
void boundary_D_1(spinor * const r, spinor * const s, su3 *u);
void boundary_D_2(spinor * const r, spinor * const s, su3 *u);
void boundary_D_3(spinor * const r, spinor * const s, su3 *u);
void boundary_D_4(spinor * const r, spinor * const s, su3 *u);
void boundary_D_5(spinor * const r, spinor * const s, su3 *u);
void boundary_D_6(spinor * const r, spinor * const s, su3 *u);
void boundary_D_7(spinor * const r, spinor * const s, su3 *u);

#endif
