/**********************************************************************
 *
 * Copyright (C) 2013  Florian Burger
 *
 * A 32-bit version of the Half-spinor implementation by Carsten Urbach
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

#ifndef _HALFSPINOR_HOPPING32_H
#define _HALFSPINOR_HOPPING32_H

#ifdef _prefetch_spinor
#undef _prefetch_spinor
#endif
#define _prefetch_spinor(s)
#ifdef _prefetch_halfspinor
#undef _prefetch_halfspinor
#endif
#define _prefetch_halfspinor(hs)
#ifdef _prefetch_spinor_32
#undef _prefetch_spinor_32
#endif
#define _prefetch_spinor_32(s)
#ifdef _prefetch_su3_32
#undef _prefetch_su3_32
#endif
#define _prefetch_su3_32(U)

#define _hop_t_p_pre32()                            \
  _vector_assign(rs.s0, s->s0);                     \
  _vector_assign(rs.s1, s->s1);                     \
  _vector_assign(rs.s2, s->s2);                     \
  _vector_assign(rs.s3, s->s3);                     \
  _vector_add(psi, rs.s0, rs.s2);                   \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s0, ka0_32, chi); \
  _vector_add(psi, rs.s1, rs.s3);                   \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s1, ka0_32, chi);

#define _hop_t_m_pre32()                   \
  _vector_sub(phi2[ix]->s0, rs.s0, rs.s2); \
  _vector_sub(phi2[ix]->s1, rs.s1, rs.s3);

#define _hop_x_p_pre32()                            \
  _vector_i_add(psi, rs.s0, rs.s3);                 \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s0, ka1_32, chi); \
  _vector_i_add(psi, rs.s1, rs.s2);                 \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s1, ka1_32, chi);

#define _hop_x_m_pre32()                     \
  _vector_i_sub(phi2[ix]->s0, rs.s0, rs.s3); \
  _vector_i_sub(phi2[ix]->s1, rs.s1, rs.s2);

#define _hop_y_p_pre32()                            \
  _vector_add(psi, rs.s0, rs.s3);                   \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s0, ka2_32, chi); \
  _vector_sub(psi, rs.s1, rs.s2);                   \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s1, ka2_32, chi);

#define _hop_y_m_pre32()                   \
  _vector_sub(phi2[ix]->s0, rs.s0, rs.s3); \
  _vector_add(phi2[ix]->s1, rs.s1, rs.s2);

#define _hop_z_p_pre32()                            \
  _vector_i_add(psi, rs.s0, rs.s2);                 \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s0, ka3_32, chi); \
  _vector_i_sub(psi, rs.s1, rs.s3);                 \
  _su3_multiply(chi, (*U), psi);                    \
  _complex_times_vector(phi2[ix]->s1, ka3_32, chi);

#define _hop_z_m_pre32()                     \
  _vector_i_sub(phi2[ix]->s0, rs.s0, rs.s2); \
  _vector_i_add(phi2[ix]->s1, rs.s1, rs.s3);

#define _hop_t_p_post32()              \
  ;                                    \
  _vector_assign(rs.s0, phi2[ix]->s0); \
  _vector_assign(rs.s2, phi2[ix]->s0); \
  _vector_assign(rs.s1, phi2[ix]->s1); \
  _vector_assign(rs.s3, phi2[ix]->s1);

#define _hop_t_m_post32()                     \
  ;                                           \
  _vector_assign(psi, phi2[ix]->s0);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka0_32, chi); \
  _vector_add_assign(rs.s0, psi);             \
  _vector_sub_assign(rs.s2, psi);             \
  _vector_assign(psi, phi2[ix]->s1);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka0_32, chi); \
  _vector_add_assign(rs.s1, psi);             \
  _vector_sub_assign(rs.s3, psi);

#define _hop_x_p_post32()                    \
  ;                                          \
  _vector_add_assign(rs.s0, phi2[ix]->s0);   \
  _vector_i_sub_assign(rs.s3, phi2[ix]->s0); \
  _vector_add_assign(rs.s1, phi2[ix]->s1);   \
  _vector_i_sub_assign(rs.s2, phi2[ix]->s1);

#define _hop_x_m_post32()                     \
  ;                                           \
  _vector_assign(psi, phi2[ix]->s0);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka1_32, chi); \
  _vector_add_assign(rs.s0, psi);             \
  _vector_i_add_assign(rs.s3, psi);           \
  _vector_assign(psi, phi2[ix]->s1);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka1_32, chi); \
  _vector_add_assign(rs.s1, psi);             \
  _vector_i_add_assign(rs.s2, psi);

#define _hop_y_p_post32()                  \
  ;                                        \
  _vector_add_assign(rs.s0, phi2[ix]->s0); \
  _vector_add_assign(rs.s3, phi2[ix]->s0); \
  _vector_add_assign(rs.s1, phi2[ix]->s1); \
  _vector_sub_assign(rs.s2, phi2[ix]->s1);

#define _hop_y_m_post32()                     \
  ;                                           \
  _vector_assign(psi, phi2[ix]->s0);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka2_32, chi); \
  _vector_add_assign(rs.s0, psi);             \
  _vector_sub_assign(rs.s3, psi);             \
  _vector_assign(psi, phi2[ix]->s1);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka2_32, chi); \
  _vector_add_assign(rs.s1, psi);             \
  _vector_add_assign(rs.s2, psi);

#define _hop_z_p_post32()                    \
  ;                                          \
  _vector_add_assign(rs.s0, phi2[ix]->s0);   \
  _vector_i_sub_assign(rs.s2, phi2[ix]->s0); \
  _vector_add_assign(rs.s1, phi2[ix]->s1);   \
  _vector_i_add_assign(rs.s3, phi2[ix]->s1);

#define _hop_z_m_post32()                     \
  ;                                           \
  _vector_assign(psi, phi2[ix]->s0);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka3_32, chi); \
  _vector_add_assign(rs.s0, psi);             \
  _vector_i_add_assign(rs.s2, psi);           \
  _vector_assign(psi, phi2[ix]->s1);          \
  _su3_inverse_multiply(chi, (*U), psi);      \
  _complexcjg_times_vector(psi, ka3_32, chi); \
  _vector_add_assign(rs.s1, psi);             \
  _vector_i_sub_assign(rs.s3, psi);

#define _hop_mul_g5_cmplx_and_store32(res)             \
  _complex_times_vector((res)->s0, cfactor, rs.s0);    \
  _complex_times_vector((res)->s1, cfactor, rs.s1);    \
  _complexcjg_times_vector((res)->s2, cfactor, rs.s2); \
  _complexcjg_times_vector((res)->s3, cfactor, rs.s3);

#define _g5_cmplx_sub_hop_and_g5store32(res)       \
  _complex_times_vector(psi, cfactor, pn->s0);     \
  _vector_sub((res)->s0, psi, rs.s0);              \
  _complex_times_vector(psi2, cfactor, pn->s1);    \
  _vector_sub((res)->s1, psi2, rs.s1);             \
  _complexcjg_times_vector(psi, cfactor, pn->s2);  \
  _vector_sub((res)->s2, rs.s2, psi);              \
  _complexcjg_times_vector(psi2, cfactor, pn->s3); \
  _vector_sub((res)->s3, rs.s3, psi2);

#define _hop_store_post32(res)    \
  _vector_assign(res->s0, rs.s0); \
  _vector_assign(res->s1, rs.s1); \
  _vector_assign(res->s2, rs.s2); \
  _vector_assign(res->s3, rs.s3);

#define _declare_hregs() \
  spinor32 ALIGN32 rs;   \
  su3_vector32 ALIGN32 psi, chi;

//, psi2, chi2;

#endif
