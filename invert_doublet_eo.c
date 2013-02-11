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
 *
 * invert_doublet_eo makes an inversion with EO precoditioned
 * tm Operator with a nondegenerate doublet
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_doublet_eo returns the number of iterations neede or -1 if the 
 * solver did not converge.
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ****************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include<stdlib.h>
#include"global.h"
#include"linalg_eo.h"
#include"operator/tm_operators.h"
#include"operator/Hopping_Matrix.h"
#include"operator/D_psi.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange/xchange.h"
#include"operator/tm_operators_nd.h"
#include"invert_doublet_eo.h"


#ifdef HAVE_GPU
#  include"GPU/cudadefs.h"
#  include"temporalgauge.h"
#  include"measure_gauge_action.h"
int mixedsolve_eo_nd (spinor *, spinor *, spinor *, spinor *, int, double, int);
int mixedsolve_eo_nd_mpi(spinor *, spinor *, spinor *, spinor *, int, double, int);
#  ifdef TEMPORALGAUGE
extern su3* g_trafo;
#  endif
#endif


int invert_doublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
		      spinor * const Even_new_c, spinor * const Odd_new_c, 
		      spinor * const Even_s, spinor * const Odd_s,
		      spinor * const Even_c, spinor * const Odd_c,
		      const double precision, const int max_iter,
		      const int solver_flag, const int rel_prec) {

  int iter = 0;
  
  
#ifdef HAVE_GPU
#  ifdef TEMPORALGAUGE
  
  /* initialize temporal gauge here */
  int retval;
  double dret1, dret2;
  double plaquette1 = 0.0;
  double plaquette2 = 0.0;
  
  if (usegpu_flag) {
    
    /* need VOLUME here (not N=VOLUME/2)*/
    if ((retval = init_temporalgauge_trafo(VOLUME, g_gauge_field)) != 0 ) {				// initializes the transformation matrices
      if (g_proc_id == 0) printf("Error while gauge fixing to temporal gauge. Aborting...\n");   	//	g_tempgauge_field as a copy of g_gauge_field
      exit(200);
    }
    
    /* do trafo */
    plaquette1 = measure_plaquette(g_gauge_field);
    apply_gtrafo(g_gauge_field, g_trafo);								// transformation of the gauge field
    plaquette2 = measure_plaquette(g_gauge_field);
    if (g_proc_id == 0) printf("\tPlaquette before gauge fixing: %.16e\n", plaquette1/6./VOLUME);
    if (g_proc_id == 0) printf("\tPlaquette after gauge fixing:  %.16e\n", plaquette2/6./VOLUME);
    
    /* do trafo to odd_s part of source */
    dret1 = square_norm(Odd_s, VOLUME/2 , 1);
    apply_gtrafo_spinor_odd(Odd_s, g_trafo);								// odd spinor transformation, strange
    dret2 = square_norm(Odd_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    /* do trafo to odd_c part of source */
    dret1 = square_norm(Odd_c, VOLUME/2 , 1);
    apply_gtrafo_spinor_odd(Odd_c, g_trafo);								// odd spinor transformation, charm
    dret2 = square_norm(Odd_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);       
    
    /* do trafo to even_s part of source */
    dret1 = square_norm(Even_s, VOLUME/2 , 1);
    apply_gtrafo_spinor_even(Even_s, g_trafo);							// even spinor transformation, strange
    dret2 = square_norm(Even_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    /* do trafo to even_c part of source */
    dret1 = square_norm(Even_c, VOLUME/2 , 1);
    apply_gtrafo_spinor_even(Even_c, g_trafo);							// even spinor transformation, charm
    dret2 = square_norm(Even_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
#    ifdef MPI
    xchange_gauge(g_gauge_field);
#    endif
    
  } 
#  endif  
#endif /* HAVE_GPU*/


  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using even/odd preconditioning!\n"); fflush(stdout);}
  M_ee_inv_ndpsi(Even_new_s, Even_new_c, 
		 Even_s, Even_c,
		 g_mubar, g_epsbar);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new_s);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI+1], Even_new_c);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd_s, VOLUME/2);
  assign_mul_add_r(g_spinor_field[DUM_DERI+1], +1., Odd_c, VOLUME/2);
  
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  
  /* Here we invert the hermitean operator squared */
  
  if(g_proc_id == 0) {
    printf("# Using CG for TMWILSON flavour doublet!\n"); 
    fflush(stdout);
  }
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
  
  
#ifdef HAVE_GPU
  if (usegpu_flag) {	// GPU, mixed precision solver
#  if defined(MPI) && defined(PARALLELT)
    iter = mixedsolve_eo_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
			    max_iter, precision, rel_prec);
#  elif !defined(MPI) && !defined(PARALLELT)
    iter = mixedsolve_eo_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
			    max_iter, precision, rel_prec);
#  else
    printf("MPI and/or PARALLELT are not appropriately set for the GPU implementation. Aborting...\n");
    exit(-1);
#  endif
  }
  else {		// CPU, conjugate gradient
    iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		     max_iter, precision, rel_prec, 
		     VOLUME/2, &Qtm_pm_ndpsi);
  }
#else			// CPU, conjugate gradient
  iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		   max_iter, precision, rel_prec, 
		   VOLUME/2, &Qtm_pm_ndpsi);
#endif
  
  
  Qtm_dagger_ndpsi(Odd_new_s, Odd_new_c,
		   Odd_new_s, Odd_new_c);

  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  M_ee_inv_ndpsi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
		 g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		 g_mubar, g_epsbar);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
  
  
#ifdef HAVE_GPU  
  /* return from temporal gauge again */
#  ifdef TEMPORALGAUGE
  
  if (usegpu_flag) { 
    
    /* undo trafo */
    /* apply_inv_gtrafo(g_gauge_field, g_trafo);*/
    /* copy back the saved original field located in g_tempgauge_field -> update necessary*/
    plaquette1 = measure_plaquette(g_gauge_field);
    copy_gauge_field(g_gauge_field, g_tempgauge_field);
    g_update_gauge_copy = 1;
    plaquette2 = measure_plaquette(g_gauge_field);
    if (g_proc_id == 0) printf("\tPlaquette before inverse gauge fixing: %.16e\n", plaquette1/6./VOLUME);
    if (g_proc_id == 0) printf("\tPlaquette after inverse gauge fixing:  %.16e\n", plaquette2/6./VOLUME);
    
    /* undo trafo to source Even_s */
    dret1 = square_norm(Even_s, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_even(Even_s, g_trafo);
    dret2 = square_norm(Even_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    
    /* undo trafo to source Even_c */
    dret1 = square_norm(Even_c, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_even(Even_c, g_trafo);
    dret2 = square_norm(Even_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1);
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2); 
    
    /* undo trafo to source Odd_s */
    dret1 = square_norm(Odd_s, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_odd(Odd_s, g_trafo);
    dret2 = square_norm(Odd_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    /* undo trafo to source Odd_c */
    dret1 = square_norm(Odd_c, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_odd(Odd_c, g_trafo);
    dret2 = square_norm(Odd_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2); 
    
    
    // Even_new_s
    dret1 = square_norm(Even_new_s, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_even(Even_new_s, g_trafo);
    dret2 = square_norm(Even_new_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    // Even_new_c
    dret1 = square_norm(Even_new_c, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_even(Even_new_c, g_trafo);
    dret2 = square_norm(Even_new_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    // Odd_new_s
    dret1 = square_norm(Odd_new_s, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_odd(Odd_new_s, g_trafo);
    dret2 = square_norm(Odd_new_s, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2);
    
    // Odd_new_c
    dret1 = square_norm(Odd_new_c, VOLUME/2 , 1);
    apply_inv_gtrafo_spinor_odd(Odd_new_c, g_trafo);
    dret2 = square_norm(Odd_new_c, VOLUME/2, 1);
    if (g_proc_id == 0) printf("\tsquare norm before gauge fixing: %.16e\n", dret1); 
    if (g_proc_id == 0) printf("\tsquare norm after gauge fixing:  %.16e\n", dret2); 
    
    finalize_temporalgauge();
    
#    ifdef MPI
    xchange_gauge(g_gauge_field);
#    endif
    
  }
#  endif
#endif
  return(iter);
}


int invert_cloverdoublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
			    spinor * const Even_new_c, spinor * const Odd_new_c, 
			    spinor * const Even_s, spinor * const Odd_s,
			    spinor * const Even_c, spinor * const Odd_c,
			    const double precision, const int max_iter,
			    const int solver_flag, const int rel_prec) {
  
  int iter = 0;
  
  
  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using even/odd preconditioning!\n"); fflush(stdout);}
  Msw_ee_inv_ndpsi(Even_new_s, Even_new_c, 
		   Even_s, Even_c);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI], Even_new_s);
  Hopping_Matrix(OE, g_spinor_field[DUM_DERI+1], Even_new_c);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_mul_add_r(g_spinor_field[DUM_DERI], +1., Odd_s, VOLUME/2);
  assign_mul_add_r(g_spinor_field[DUM_DERI+1], +1., Odd_c, VOLUME/2);
  
  /* Do the inversion with the preconditioned  */
  /* matrix to get the odd sites               */
  
  /* Here we invert the hermitean operator squared */
  
  if(g_proc_id == 0) {
    printf("# Using CG for TMWILSON flavour doublet!\n"); 
    fflush(stdout);
  }
  gamma5(g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI], VOLUME/2);
  gamma5(g_spinor_field[DUM_DERI+1], g_spinor_field[DUM_DERI+1], VOLUME/2);
  
  iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		   max_iter, precision, rel_prec, 
		   VOLUME/2, &Qsw_pm_ndpsi);
  
  
  Qsw_dagger_ndpsi(Odd_new_s, Odd_new_c,
		   Odd_new_s, Odd_new_c);
  
  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  Msw_ee_inv_ndpsi(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
		   g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);
  
  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
  
  return(iter);
}

