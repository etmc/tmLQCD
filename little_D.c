/* $Id$ */
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "complex.h"
#include "deflation/deflation_block.h"
#include "linalg/blas.h"
#include "little_D.h"


/* assume we have a little field w                       */
/* which has length 9*no_blocks*N_s                      */
/* with usual order in space                             */
/* no_blocks = 2 currently fixed                         */
/* and blocks devide z-direction by 2                    */
/*                                                       */
/* block[0], block[1], block[0], block[1], block[0]  ... */
/* local             , +t                , -t        ... */

const int no_blocks = 2;
const int nblks_t = 1;
const int nblks_x = 1;
const int nblks_y = 1;
const int nblks_z = 2;
int nblks_dir[4] = {1,1,1,2};

/* some lapack related stuff */
static int ONE = 1;
static complex CONE, CZERO, CMONE;

void init_little_field_exchange(complex * w);
void wait_little_field_exchange(const int i, complex * w);


void little_D(complex * v, complex *w) {
  int i, j, k, sq = g_N_s*g_N_s;
  CONE.re = 1.;
  CONE.im = 0.;
  CMONE.re = -1.;
  CMONE.im = 0.;
  CZERO.re = 0.;
  CZERO.im = 0.;
  
#ifdef MPI
  init_little_field_exchange(w);
#endif
  /* all the mpilocal stuff first */
  for(i = 0; i < no_blocks; i++) {
    /* diagonal term */
    _FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, g_blocks[i].little_dirac_operator,
	       &g_N_s, w + i*g_N_s, &ONE, &CZERO, v + i*g_N_s, &ONE, 1);

    /* offdiagonal terms */
    for(j = 0; j < 8; j++) {
      /* set k to neighbour in direction j */
      k = g_blocks[i].mpilocal_neighbour[j];
      /* if k is on the same mpi proc, but not myself */
      if(k > -1 && k != i) {
	_FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, g_blocks[i].little_dirac_operator + (j+1)*sq,
		   &g_N_s, w + k*g_N_s, &ONE, &CONE, v + i*g_N_s, &ONE, 1);
      }
    }
  }

#ifdef MPI
  /* now all non-mpilocal stuff */
  /* start with z direction     */
  for(j = 8; j > -1; j--) {
    wait_little_field_exchange(j, w);
    for(i = 0; i < no_blocks; i++) {
      k = g_blocks[i].mpilocal_neighbour[j];
      if(k < 0) {
	_FT(zgemv)("N", &g_N_s, &g_N_s, &CONE, g_blocks[i].little_dirac_operator + (j+1)*sq,
		   &g_N_s, w + ((j+1)*no_blocks + i)*g_N_s, &ONE, &CONE, v + i*g_N_s, &ONE, 1);
      }
    }
  }
  
#endif
  return;
}
