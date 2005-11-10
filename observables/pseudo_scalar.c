/**************************************************************************
 *
 * The externally accessible function is
 *
 *   void calculate_pseudo_scalar_prop(spinor *** const qprop, int nstore, int im)
 *
 *     Calculates the 2-point function of the pseudo-scalar density
 *
 *                 psi_bar gamma_5 psi
 *
 *     using the quark propagator qprop for mass im
 *
 * Author: Ines Wetzorke <Ines.Wetzorke@desy.de> Feb 2003
 *
 *
 *   The momentum of the meson can now be also nonzero
 *
 *                Stefano Capitani <stefano@ifh.de>, Apr 2004
 *
 **************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "bin/global.h"
#include "bin/io.h"
#include "sse-extension/sse.h"
#include "su3/su3.h"
#include "su3/su3adj.h"
#include "geometry_eo.h"
#include "pseudo_scalar.h"

void pseudo_scalar(spinor * const qprop, int nstore, int im){

  int x0,x1,x2,x3,ix;
  int global_x1, global_x2, global_x3;
  double pp, pprop[T];

  FILE *ofs;
  char *filename_stub = "pssca_corr_";
  char *filename;
  char buf[100];

  filename=buf;
  sprintf(filename,"%s%d.dat", filename_stub,im);

  for (x0=0;x0<T;x0++){
    ks=0.0;
    kc=0.0;
    for (x1=0;x1<LX;x1++){
      for (x2=0;x2<LY;x2++){
	for (x3=0;x3<LZ;x3++){
	  ix=index(x0,x1,x2,x3);
          global_x1=x1+g_proc_coords[0]*LX;
          global_x2=x2+g_proc_coords[1]*LY;
          global_x3=x3+g_proc_coords[2]*LZ;
	  for(is=0;is<4;is++){
	    for(ic=0;ic<3;ic++){
	      pp=_spinor_prod_re(qprop[is][ic][ix],qprop[is][ic][ix]);
              pp=pp*cos(g_p1*global_x1+g_p2*global_x2+g_p3*global_x3);
	      tr=pp+kc;
	      ts=tr+ks;
	      tt=ts-ks;
	      ks=ts;
	      kc=tr-tt;
	    } /* ic */
	  } /* is */
	} /* x3 */
      } /* x2 */
    } /* x1 */
#if defined MPI
    kc=ks+kc;
    MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, g_cart_grid);
    pprop[x0]=ks;
#else
    pprop[x0]=ks+kc;
#endif
  } /* x0 */

#if defined MPI
  if (g_proc_id == 0) {
#endif
  if (nstore == 0) {
    ofs = fopen(filename,"w");
    write_info(ofs,im);
    fprintf(ofs, "# T  pseudo-scalar\n");
  }
  else {
    ofs = fopen(filename,"a");
  }
  fprintf(ofs, "# measurement %d\n", nstore+1);
  for (x0=0;x0<T;x0++){
    fprintf(ofs, "%3d  %+12.10e\n", x0, pprop[x0]); 
  }
  fclose(ofs);
#if defined MPI
  }
#endif
}
