/* $Id$ */
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"global.h"
#include"su3.h"
#include"complex.h"
#include"start.h"
#include"D_psi.h"
#include"poly_precon.h"
#include"linalg_eo.h"
#include"gram-schmidt.h"
#include"generate_dfl_subspace.h"

spinor ** dfl_fields;


int generate_dfl_subspace(const int Ns, const int N) {
  int i,j, vpr = VOLUMEPLUSRAND*sizeof(spinor)/sizeof(complex), 
    vol = VOLUME*sizeof(spinor)/sizeof(complex);
  double nrm, e = 0.5, d = 1.;
  
  for(i = 0; i < Ns; i++) {
    random_spinor_field(dfl_fields[i], 1, N);

    ModifiedGS((complex*)dfl_fields[i], vol, i, (complex*)dfl_fields[0], vpr);
    nrm = sqrt(square_norm(dfl_fields[i], N));
    mul_r(dfl_fields[i], 1./nrm, dfl_fields[i], N);

    for(j = 0; j < 10; j++) {
      poly_nonherm_precon(g_spinor_field[DUM_SOLVER], dfl_fields[i], e, d, 20, N);

      ModifiedGS((complex*)g_spinor_field[DUM_SOLVER], vol, i, (complex*)dfl_fields[0], vpr);
      nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N));
      mul_r(dfl_fields[i], 1./nrm, g_spinor_field[DUM_SOLVER], N);
    }
    /* test quality */
    D_psi(g_spinor_field[DUM_SOLVER], dfl_fields[i]);
    nrm = sqrt(square_norm(g_spinor_field[DUM_SOLVER], N));
    if(g_proc_id == 0 && g_debug_level > -1) {
      printf(" ||D psi_%d||/||psi_%d|| = %1.5e\n", i, i, nrm); 
    }
  }
  
  return(0);
}
