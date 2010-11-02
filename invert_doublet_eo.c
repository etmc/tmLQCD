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
/* $Id$ */

/****************************************************************
 *
 * invert_eo makes an inversion with EO precoditioned
 * tm Operator
 *
 * Even and Odd are the numbers of spinor_field that contain
 * the even and the odd sites of the source. The result is stored
 * int Even_new and Odd_new.
 *
 * invert_eo returns the number of iterations neede or -1 if the 
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
#include"tm_operators.h"
#include"Hopping_Matrix.h"
#include"D_psi.h"
#include"linsolve.h"
#include"gamma.h"
#include"solver/solver.h"
#include"read_input.h"
#include"xchange.h"
#include"Nondegenerate_Matrix.h"
#include"invert_doublet_eo.h"


#ifdef HAVE_GPU
  #include"GPU/cudadefs.h"
  #include"temporalgauge.h"
  #include"observables.h"
  int mixedsolve_eo_nd (spinor *, spinor *, spinor *, spinor *, int, double, int);
  int mixedsolve_eo_nd_mpi(spinor *, spinor *, spinor *, spinor *, int, double, int);
  #ifdef TEMPORALGAUGE
    extern su3* g_trafo;
  #endif
#endif


int invert_doublet_eo(spinor * const Even_new_s, spinor * const Odd_new_s, 
		      spinor * const Even_new_c, spinor * const Odd_new_c, 
		      spinor * const Even_s, spinor * const Odd_s,
		      spinor * const Even_c, spinor * const Odd_c,
		      const double precision, const int max_iter,
		      const int solver_flag, const int rel_prec) {

  int iter = 0;
  
  
  #ifdef HAVE_GPU
  #ifdef TEMPORALGAUGE
    /* initialize temporal gauge here */
    int retval;
    double dret;
    double plaquette = 0.0;

    if(usegpu_flag){
    
      /* need VOLUME here (not N=VOLUME/2)*/
      if((retval=init_temporalgauge_trafo(VOLUME, g_gauge_field)) !=0){
	if(g_proc_id == 0) printf("Error while gauge fixing to temporal gauge. Aborting...\n");   
	exit(200);
      }
      plaquette = measure_gauge_action();
      if(g_proc_id == 0) printf("Plaquette before gauge fixing: %.16e\n", plaquette/6./VOLUME);
      /* do trafo */
      apply_gtrafo(g_gauge_field, g_trafo);
      plaquette = measure_gauge_action();
      if(g_proc_id == 0) printf("Plaquette after gauge fixing: %.16e\n", plaquette/6./VOLUME);
    
      /* do trafo to odd_s part of source */
      dret = square_norm(Odd_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_odd(Odd_s, g_trafo);
      dret = square_norm(Odd_s, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);
      /* do trafo to odd_c part of source */
      dret = square_norm(Odd_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_odd(Odd_c, g_trafo);
      dret = square_norm(Odd_c, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);       
    
      /* do trafo to even_s part of source */
      dret = square_norm(Even_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_even(Even_s, g_trafo);
      dret = square_norm(Even_s, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret);
      /* do trafo to even_c part of source */
      dret = square_norm(Even_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_gtrafo_spinor_even(Even_c, g_trafo);
      dret = square_norm(Even_c, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret);      
    } 
#endif  
#endif /* HAVE_GPU*/


  /* here comes the inversion using even/odd preconditioning */
  if(g_proc_id == 0) {printf("# Using Even/Odd preconditioning!\n"); fflush(stdout);}
  M_ee_inv_ND(Even_new_s, Even_new_c, 
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
  
  
  #ifdef HAVE_GPU
    if (usegpu_flag) {	// GPU, mixed precision solver
      #if defined(MPI) && defined(PARALLELT)
        iter = mixedsolve_eo_nd_mpi(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
                                    max_iter, precision, rel_prec);
      #elif !defined(MPI) && !defined(PARALLELT)
        iter = mixedsolve_eo_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
                                max_iter, precision, rel_prec);
      #else
        printf("MPI and/or PARALLELT are not appropriately set for the GPU implementation. Aborting...\n");
        exit(-1);
      #endif
    }
    else {		// CPU, conjugate gradient
      iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		       max_iter, precision, rel_prec, 
		       VOLUME/2, &Q_Qdagger_ND, 0, 1000);
    }
  #else			// CPU, conjugate gradient
    iter = cg_her_nd(Odd_new_s, Odd_new_c, g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1],
		     max_iter, precision, rel_prec, 
		     VOLUME/2, &Q_Qdagger_ND, 0, 1000);
  #endif
  
  
  QdaggerNon_degenerate(Odd_new_s, Odd_new_c,
			Odd_new_s, Odd_new_c);
  
  /* Reconstruct the even sites                */
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI], Odd_new_s);
  Hopping_Matrix(EO, g_spinor_field[DUM_DERI+1], Odd_new_c);
  M_ee_inv_ND(g_spinor_field[DUM_DERI+2], g_spinor_field[DUM_DERI+3],
	      g_spinor_field[DUM_DERI], g_spinor_field[DUM_DERI+1]);

  /* The sign is plus, since in Hopping_Matrix */
  /* the minus is missing                      */
  assign_add_mul_r(Even_new_s, g_spinor_field[DUM_DERI+2], +1., VOLUME/2);
  assign_add_mul_r(Even_new_c, g_spinor_field[DUM_DERI+3], +1., VOLUME/2);
  
  
  #ifdef HAVE_GPU  
    /* return from temporal gauge again */
  #ifdef TEMPORALGAUGE
    if(usegpu_flag){ 
      plaquette = measure_gauge_action();
      if(g_proc_id == 0) printf("Plaquette before inverse gauge fixing: %.16e\n", plaquette/6./VOLUME);
    
      /* undo trafo */
    
      /*apply_inv_gtrafo(g_gauge_field, g_trafo);*/
      /* copy back the saved original field located in g_tempgauge_field -> update necessary*/
      copy_gauge_field(g_gauge_field, g_tempgauge_field);
      g_update_gauge_copy = 1;
    
    
      plaquette = measure_gauge_action();
      if(g_proc_id == 0) printf("Plaquette after inverse gauge fixing: %.16e\n", plaquette/6./VOLUME);
   
      /* undo trafo to source Even_s */
      dret = square_norm(Even_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16f\n", dret); 
      apply_inv_gtrafo_spinor_even(Even_s, g_trafo);
      dret = square_norm(Even_s, VOLUME/2, 1);
      /* undo trafo to source Even_c */
      dret = square_norm(Even_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16f\n", dret); 
      apply_inv_gtrafo_spinor_even(Even_c, g_trafo);
      dret = square_norm(Even_c, VOLUME/2, 1);
      
      /* undo trafo to source Odd_s */
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret);  
      dret = square_norm(Odd_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16f\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd_s, g_trafo);
      dret = square_norm(Odd_s, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret);
      /* undo trafo to source Odd_c */
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret);  
      dret = square_norm(Odd_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16f\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd_c, g_trafo);
      dret = square_norm(Odd_c, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16f\n", dret); 
    
      // Even_new_s
      dret = square_norm(Even_new_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_even(Even_new_s, g_trafo);
      dret = square_norm(Even_new_s, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);
      // Even_new_c
      dret = square_norm(Even_new_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_even(Even_new_c, g_trafo);
      dret = square_norm(Even_new_c, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);
      
      // Odd_new_s
      dret = square_norm(Odd_new_s, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd_new_s, g_trafo);
      dret = square_norm(Odd_new_s, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret);
      // Odd_new_c
      dret = square_norm(Odd_new_c, VOLUME/2 , 1);
      if(g_proc_id == 0) printf("square norm before gauge fixing: %.16e\n", dret); 
      apply_inv_gtrafo_spinor_odd(Odd_new_c, g_trafo);
      dret = square_norm(Odd_new_c, VOLUME/2, 1);
      if(g_proc_id == 0) printf("square norm after gauge fixing: %.16e\n", dret); 
  
    
      finalize_temporalgauge();
    }
  #endif
  #endif


  return(iter);
}

