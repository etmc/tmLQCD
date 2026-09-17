
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

#include "Dov_proj.h"
#include <stdlib.h>
#include "Dov_psi.h"
#include "gamma.h"
#include "global.h"
#include "linalg_eo.h"
#include "su3.h"

void Dov_proj_plus(spinor *const R, spinor *const S) {
  spinor *aux;
  int N = VOLUMEPLUSRAND;

  aux = calloc(VOLUMEPLUSRAND, sizeof(spinor));

  Proj(aux, S, N, _PLUS);
  Dov_psi(R, aux);
  Proj(R, R, N, _PLUS);

  free(aux);
}

void Dov_proj_minus(spinor *const R, spinor *const S) {
  spinor *aux;
  int N = VOLUMEPLUSRAND;

  aux = calloc(VOLUMEPLUSRAND, sizeof(spinor));

  Proj(aux, S, N, _MINUS);
  Dov_psi(R, aux);
  Proj(R, R, N, _MINUS);

  free(aux);
}
