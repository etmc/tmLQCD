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

//This works with tm and tm+clover 
void D_psi(spinor * const P, spinor * const Q);
void D_psi_32(spinor32 * const P, spinor32 * const Q);
void D_psi_prec(spinor * const P, spinor * const Q);

//works for tm and tm+clover
void Block_D_psi(block * blk, spinor * const rr, spinor * const s);
void Block_H_psi(block * blk, spinor * const rr, spinor * const s, const int eo);

void Block_D_psi_32(block * blk, spinor32 * const rr, spinor32 * const s);
void Block_H_psi_32(block * blk, spinor32 * const rr, spinor32 * const s, const int eo);



/* csw=0 version*/
void Dtm_psi(spinor * const P, spinor * const Q);
/* csw>0 version*/
void Dsw_psi(spinor * const P, spinor * const Q);

//c_sw=0
void Block_Dtm_psi(block * blk, spinor * const rr, spinor * const s);
//c_sw > 0
void Block_Dsw_psi(block * blk, spinor * const rr, spinor * const s);

//c_sw=0
void Block_Dtm_psi_32(block * blk, spinor32 * const rr, spinor32 * const s);
//c_sw > 0
void Block_Dsw_psi_32(block * blk, spinor32 * const rr, spinor32 * const s);


void boundary_D_0(spinor * const r, spinor * const s, su3 *u);
void boundary_D_1(spinor * const r, spinor * const s, su3 *u);
void boundary_D_2(spinor * const r, spinor * const s, su3 *u);
void boundary_D_3(spinor * const r, spinor * const s, su3 *u);
void boundary_D_4(spinor * const r, spinor * const s, su3 *u);
void boundary_D_5(spinor * const r, spinor * const s, su3 *u);
void boundary_D_6(spinor * const r, spinor * const s, su3 *u);
void boundary_D_7(spinor * const r, spinor * const s, su3 *u);

#endif
