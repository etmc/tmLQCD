/* $Id$ */
#ifndef _SOURCE_GENERATION_H
#define _SOURCE_GENERATION_H

void source_generation_pion_only(spinor * const P, spinor * const Q,
				 const int t,
				 const int sample, const int nstore);

void source_generation_nucleon(spinor * const P, spinor * const Q, 
			       const int is, const int ic,
			       const int t, const int nt, const int nx, 
			       const int sample, const int nstore,
			       const int meson);

void extended_pion_source(spinor * const P, spinor * const Q,
			  spinor * const R, spinor * const S,
			  const int t0,
			  const double px, const double py, const double pz);
#endif
