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
#ifndef _SOURCE_GENERATION_H
#define _SOURCE_GENERATION_H

void gaussian_volume_source(spinor * const P, spinor * const Q,
			    const int sample, const int nstore, const int f);

void source_generation_pion_only(spinor * const P, spinor * const Q,
				 const int t, const int sample, 
                                 const int nstore, const unsigned int _seed);

void source_generation_nucleon(spinor * const P, spinor * const Q, 
			       const int is, const int ic,
			       const int t, const int nt, const int nx, 
			       const int sample, const int nstore,
			       const int meson);

void extended_pion_source(spinor * const P, spinor * const Q,
			  spinor * const R, spinor * const S,
			  const int t0, const int ts,
			  const double px, const double py, const double pz);

void source_generation_pion_zdir(spinor * const P, spinor * const Q,
                                 const int z,
                                 const int sample, const int nstore);

void full_source_spinor_field_spin_diluted_oet_ts(spinor * const full_spinor,
                                                  const int src_ts,
                                                  const int src_d,
                                                  const int sample,
                                                  const int nstore,
                                                  const unsigned int oet_seed);

void eo_source_spinor_fied_spin_diluted_oet_ts(spinor * const even_cb_spinor,
                                               spinor * const odd_cb_spinor,
                                               const int src_ts,
                                               const int src_d,
                                               const int sample,
                                               const int nstore,
                                               const unsigned int oet_seed);

// generates a point source at the global coordinates passed via the
// four element vector global_txyz_src_pos in the ordering {t,x,y,z}
// the spin index 'is' and the colour index 'ic'
void full_source_spinor_field_point(spinor * const full_spinor,
                                    const int is, const int ic,
                                    const int * const global_txyz_src_pos);

// as full_source_spinor_field_point but with output directly to checkerboarded
// spinors
void eo_source_spinor_field_point(spinor * const even_cb_spinor,
                                  spinor * const odd_cb_spinor,
                                  const int is, const int ic,
                                  const int * const global_txyz_src_pos);

#endif
