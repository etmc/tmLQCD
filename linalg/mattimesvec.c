/* $Id$ */

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "complex.h"
#include "mattimesvec.h"

/* v = M*w                         */
/* v,w complex vectors of length N */
/* M a NxN complex matrix with     */
/* leading dimension ldM >= N      */
/* we should provide special SSE2  */
/* and BG/P versions               */

void mattimesvec(complex * const v, complex * const M, complex * const w, 
		 const int N, const int ldM) {
  int i, j;

  for(i = 0; i < N; i++) {
    _mult_assign_complex(v[i], M[i*ldM], w[0]);
    for(j = 1; j < N; j++) {
      _add_assign_complex(v[i], M[i*ldM + j], w[j]);
    }
  }
  return;
}
