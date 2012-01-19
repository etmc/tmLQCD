
/**********************************************************
 *
 * Dov_proj_plus and Dov_proj_minus
 * are the projections of Dov onto the
 * positive and negative chiral sector, respectively
 *
 * Both need one work_field!
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *         Die Sep 21 15:21:33 CEST 2004
 *
 **********************************************************/

#include <stdlib.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "Dov_proj.h"
#include "gamma.h"
#include "Dov_psi.h"


void Dov_proj_plus(spinor * const R, spinor * const S)
{
  spinor *aux_ = NULL, *aux;
  int N = VOLUMEPLUSRAND;

#if ( defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif
  
  Proj(aux, S, N, _PLUS);
  Dov_psi(R, aux);
  Proj(R, R, N, _PLUS);
  
  free(aux_);
}


void Dov_proj_minus(spinor * const R, spinor * const S)
{
  spinor *aux_ = NULL, *aux;
  int N = VOLUMEPLUSRAND;

#if ( defined SSE || defined SSE2 || defined SSE3)
  aux_=calloc(VOLUMEPLUSRAND+1, sizeof(spinor));
  aux = (spinor *)(((unsigned long int)(aux_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  aux_=calloc(VOLUMEPLUSRAND, sizeof(spinor));
  aux = aux_;
#endif
  
  Proj(aux, S, N, _MINUS);
  Dov_psi(R, aux);
  Proj(R, R, N, _MINUS);
  
  free(aux_);
}



