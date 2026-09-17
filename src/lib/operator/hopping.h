/**********************************************************************
 *
 * Copyright (C) 2012 Carsten Urbach
 *
 * BG and halfspinor versions (C) 2007, 2008 Carsten Urbach
 *
 * This file is based on an implementation of the Dirac operator
 * written by Martin Luescher, modified by Martin Hasenbusch in 2002
 * and modified and extended by Carsten Urbach from 2003-2008
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
 **********************************************************************/

#ifndef _HOPPING_H
#define _HOPPING_H

#define _declare_regs()	     \
  su3_vector ALIGN psi, chi; \
  spinor ALIGN temp;

#define _hop_t_p()				\
  _vector_add(psi,sp->s0,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka0,chi);		\
  _vector_assign(temp.s0,psi);			\
  _vector_assign(temp.s2,psi);			\
  _vector_add(psi,sp->s1,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka0,chi);		\
  _vector_assign(temp.s1,psi);			\
  _vector_assign(temp.s3,psi);
  
#define _hop_t_m()				\
  _vector_sub(psi,sm->s0,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_sub_assign(temp.s2,psi);		\
  _vector_sub(psi,sm->s1,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka0,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_sub_assign(temp.s3,psi);

#define _hop_x_p()				\
  _vector_i_add(psi,sp->s0,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka1,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_sub_assign(temp.s3,psi);		\
  _vector_i_add(psi,sp->s1,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka1,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_sub_assign(temp.s2,psi);

#define _hop_x_m()				\
  _vector_i_sub(psi,sm->s0,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_add_assign(temp.s3,psi);		\
  _vector_i_sub(psi,sm->s1,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka1,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_add_assign(temp.s2,psi);

#define _hop_y_p()				\
  _vector_add(psi,sp->s0,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka2,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_add_assign(temp.s3,psi);		\
  _vector_sub(psi,sp->s1,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka2,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_sub_assign(temp.s2,psi);

#define _hop_y_m()				\
  _vector_sub(psi,sm->s0,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(temp.s0,psi);		\
  _vector_sub_assign(temp.s3,psi);		\
  _vector_add(psi,sm->s1,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka2,chi);	\
  _vector_add_assign(temp.s1,psi);		\
  _vector_add_assign(temp.s2,psi);

#define _hop_z_p()				\
  _vector_i_add(psi,sp->s0,sp->s2);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka3,chi);		\
  _vector_add_assign(temp.s0,psi);		\
  _vector_i_sub_assign(temp.s2,psi);		\
  _vector_i_sub(psi,sp->s1,sp->s3);		\
  _su3_multiply(chi,(*up),psi);			\
  _complex_times_vector(psi,ka3,chi);		\
  _vector_add_assign(temp.s1,psi);		\
  _vector_i_add_assign(temp.s3,psi);

#define _hop_z_m()				\
  _vector_i_sub(psi,sm->s0,sm->s2);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(temp.s0, psi);		\
  _vector_i_add_assign(temp.s2, psi);		\
  _vector_i_add(psi,sm->s1,sm->s3);		\
  _su3_inverse_multiply(chi,(*um),psi);		\
  _complexcjg_times_vector(psi,ka3,chi);	\
  _vector_add_assign(temp.s1, psi);		\
  _vector_i_sub_assign(temp.s3, psi);

#define _hop_mul_g5_cmplx_and_store()			\
  _complex_times_vector(rn->s0, cfactor, temp.s0);	\
  _complex_times_vector(rn->s1, cfactor, temp.s1);	\
  _complexcjg_times_vector(rn->s2, cfactor, temp.s2);	\
  _complexcjg_times_vector(rn->s3, cfactor, temp.s3);

#define _g5_cmplx_sub_hop_and_g5store()			\
  _complex_times_vector(psi, cfactor, pn->s0);		\
  _vector_sub(rn->s0, psi, temp.s0);			\
  _complex_times_vector(chi, cfactor, pn->s1);		\
  _vector_sub(rn->s1, chi, temp.s1);			\
  _complexcjg_times_vector(psi, cfactor, pn->s2);	\
  _vector_sub(rn->s2, temp.s2, psi);			\
  _complexcjg_times_vector(chi, cfactor, pn->s3);	\
  _vector_sub(rn->s3, temp.s3, chi);

#define _store_res()				\
  _vector_assign(rn->s0, temp.s0);		\
  _vector_assign(rn->s1, temp.s1);		\
  _vector_assign(rn->s2, temp.s2);		\
  _vector_assign(rn->s3, temp.s3);

#endif
