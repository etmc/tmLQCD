/***********************************************************************
 *
 * Copyright (C) 2015 Florian Burger
 * partially based on cg_mms_tm_nd.c by Andrea Shindler and Carsten Urbach
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
 * Author: 2015 Florian Burger
 * 
 * This is a Multi-Shift reliable update single/double mixed CG solver
 * it expects that the shifts fulfil
 *
 * shift[0] < shift[1] < shift{2] < ... < shift[no_shifts-1]
 *
 * in modulus. The code will use shift[i]^2, which are all >0
 *
 * parameters:
 * shifts are given to the solver in solver_params->shifts
 * number of shifts is in solver_params->no_shifts
 * the operator to invert in solver_params->M_ndpsi
 * the 32 bit operator to invert in solver_params->M_ndpsi32
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
#include "tmlqcd_config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "gamma.h"
#include "linalg_eo.h"
#include "start.h"
#include "gettime.h"
#include "solver/solver.h"
#include "solver_field.h"
#include "cg_mms_tm_nd.h"
#include "mixed_cg_mms_tm_nd.h"


static spinor32 * x_up_qmms;
static spinor32 ** mms_x_up;
static spinor32 * x_dn_qmms;
static spinor32 ** mms_x_dn;

static spinor32 * d_up_qmms;
static spinor32 ** mms_d_up;
static spinor32 * d_dn_qmms;
static spinor32 ** mms_d_dn;


static void init_mms_tm_nd_32(const unsigned int nr, const unsigned int N);
static void free_mms_tm_nd_32();

int mixed_cg_mms_tm_nd(spinor ** const Pup, spinor ** const Pdn, 
		 spinor * const Qup, spinor * const Qdn, 
		 solver_params_t * solver_params) {

  double eps_sq = solver_params->squared_solver_prec;
  int noshifts = solver_params->no_shifts;
  int rel_prec = solver_params->rel_prec;
  int max_iter = solver_params->max_iter;
  int check_abs, check_rel;
  double * shifts = solver_params->shifts;
  int Nshift = noshifts;
 
  // algorithm
  double rr_up, rr_dn, rr, rr_old, r0r0, dAd_up, dAd_dn, dAd;  
  
  if(rel_prec){
    check_rel = 1;
    check_abs = 0;
   }
   else{
    check_rel = 0;
    check_abs = 1;     
  }
  
  int use_eo=1, eofactor=2;
  //not even-odd?
  if(solver_params->sdim == VOLUME) {
    eofactor = 1;
    use_eo = 0;
  }
  
  int N = VOLUME/eofactor;
  int Vol = VOLUMEPLUSRAND/eofactor;
 
  
  // norm of source
  rr_up = square_norm(Qup, N, 1);
  rr_dn = square_norm(Qdn, N, 1);
  rr    = rr_up + rr_dn;  
 
  if( (g_cart_id == 0 && g_debug_level > 2)) printf("# CGMMSND_mixed: Initial mms residue: %.6e\n", rr);
  if(rr < 1.0e-4){
    if( (g_cart_id == 0 && g_debug_level > 2)) printf("# CGMMSND_mixed: norm of source too low: falling back to double mms solver %.6e\n", rr);
    return(cg_mms_tm_nd(Pup, Pdn, Qup, Qdn, solver_params));
  }
  
  r0r0   = rr;	// for relative precision 
  rr_old = rr;	// for the first iteration
  
  
  
  //allocate an auxiliary solver fields 
  spinor ** sf = NULL;
  const int nr_sf = 6;
  init_solver_field(&sf, Vol, nr_sf);  
   
  spinor32 ** sf32 = NULL;
  const int nr_sf32 = 8;
  init_solver_field_32(&sf32, Vol, nr_sf32);  
  
  
  //spinor fields  
  //we need one less than shifts, since one field is cared of by the usual cg fields
  init_mms_tm_nd_32(noshifts-1, Vol);
   
  // Pup/dn  can be used as auxiliary field to work on, as it is not later used (could be used as initial guess at the very start)
  // Q_up/dn  can be used as feedback, or if not, also as auxiliary field
  

  
  //allocate cg constants
  double * sigma;
  double * zitam1, * zita;
  double * alphas, * betas;
  double gamma;
  double alpham1;
    sigma = (double*)calloc((noshifts), sizeof(double));
    zitam1 = (double*)calloc((noshifts), sizeof(double));
    zita = (double*)calloc((noshifts), sizeof(double));
    alphas = (double*)calloc((noshifts), sizeof(double));
    betas = (double*)calloc((noshifts), sizeof(double));



  spinor32 *  r_up, *  r_dn, * Ad_up, * Ad_dn, *  x_up, *  x_dn, *  d_up, *  d_dn;		
  spinor * r_up_d, * r_dn_d, * x_up_d, * x_dn_d, * Ax_up_d, * Ax_dn_d;
  
 // iteration counter
 int j; 
 
 //reliable update flag
 int rel_update = 0;
 //no of reliable updates done
 int no_rel_update = 0;
 //use reliable update flag
 int use_reliable = 1;
 
 double rel_delta = 1.0e-10;
 int trigger_shift = -1;
 double * res;
 double * res0;
 double * maxres;
 res = (double*)calloc((noshifts), sizeof(double));
 res0 = (double*)calloc((noshifts), sizeof(double));
 maxres = (double*)calloc((noshifts), sizeof(double)); 
    
  /////////////////
  // ASSIGNMENTS //
  /////////////////
  
  x_up  = sf32[0];	
  x_dn  = sf32[1];	
  r_up  = sf32[2];	
  r_dn  = sf32[3];
  d_up  = sf32[4];
  d_dn  = sf32[5];
  Ad_up = sf32[6];
  Ad_dn = sf32[7];


  x_up_d = sf[0];
  x_dn_d = sf[1];
  r_up_d = sf[2];
  r_dn_d = sf[3];
  Ax_up_d = sf[4];
  Ax_dn_d = sf[5];  
  
  /*
  //matrix test
   spinor32 * help_low_up = sf32[0];
   spinor32 * help_low_dn = sf32[1];   
   spinor * help_high_up = sf[0];
   spinor * help_high_dn = sf[1];   
   assign_to_32(help_low_up, Qup, N);
   assign_to_32(help_low_dn, Qdn, N);   
   assign(help_high_up, Qup, N);
   assign(help_high_dn, Qdn, N);   
   double sqn_high = square_norm(help_high_up,N,1) +
                     square_norm(help_high_dn,N,1);
   printf("square_norm(Q_high) = %e\n", sqn_high);
   float sqn_low  = square_norm_32(help_low_up,N,1) +
                    square_norm_32(help_low_dn,N,1);   
   printf("square_norm(Q_low) = %e\n", sqn_low);  
   
   solver_params->M_ndpsi32(sf32[2], sf32[3], help_low_up, help_low_dn);
   solver_params->M_ndpsi(sf[2], sf[3], help_high_up, help_high_dn);
   
   assign_to_64(sf[4], sf32[2], N);
   assign_to_64(sf[5], sf32[3], N);   
   diff(sf[0], sf[4], sf[2], N);
   diff(sf[1], sf[5], sf[3], N);   
   double sqnrm = square_norm(sf[0], N, 1) +
                  square_norm(sf[1], N, 1);
   printf("Operator 32 test: (square_norm) / (spinor component) = %.8e\n", sqnrm/24.0/N);
   exit(1);  
  */
  
  // r(0) = b
  assign_to_32(r_up, Qup, N);
  assign_to_32(r_dn, Qdn, N); 
  
  // d(0) = b
  assign_to_32(d_up, Qup, N);
  assign_to_32(d_dn, Qdn, N); 
  

  
  maxres[0] = rr;
  res[0] = rr;
  res0[0] = rr;
  alphas[0] = 1.0;
  betas[0] = 0.0;
  sigma[0] = shifts[0]*shifts[0];
  if(g_cart_id == 0 && g_debug_level > 2) printf("# CGMMSND_mixed: shift %d is %e\n", 0, sigma[0]);

  // currently only implemented for P=0 
  for(int im = 1; im < noshifts; im++) {
    maxres[im] = rr;
    res[im] = rr;
    res0[im] = rr;    
    sigma[im] = shifts[im]*shifts[im] - sigma[0];
    if(g_cart_id == 0 && g_debug_level > 2) printf("# CGMMSND_mixed: shift %d is %e\n", im, sigma[im]);
    // these will be the result spinor fields
    zero_spinor_field_32(mms_x_up[im-1], N);
    zero_spinor_field_32(mms_x_dn[im-1], N);    

    assign_to_32(mms_d_up[im-1], Qup, N);
    assign_to_32(mms_d_dn[im-1], Qdn, N);
    zitam1[im] = 1.0;
    zita[im] = 1.0;
    alphas[im] = 1.0;
    betas[im] = 0.0;
  }

  //zero fields for solution Pup, Pdn
  for(int im = 0; im < noshifts; im++){
    zero_spinor_field(Pup[im], N);
    zero_spinor_field(Pdn[im], N);    
  }
  
  
  //////////
  // LOOP //
  //////////
    
  for (j = 0; j < max_iter; j++) {   
      // A*d(k)
    solver_params->M_ndpsi32(Ad_up, Ad_dn, d_up,  d_dn);     
    //add zero'th shift
    assign_add_mul_r_32(Ad_up, d_up, (float) sigma[0], N);
    assign_add_mul_r_32(Ad_dn, d_dn, (float) sigma[0], N);
	     
    
    // alpha = r(k)*r(k) / d(k)*A*d(k)
    dAd_up = scalar_prod_r_32(d_up, Ad_up, N, 1);
    dAd_dn = scalar_prod_r_32(d_dn, Ad_dn, N, 1);

    dAd    = dAd_up + dAd_dn; 
    alpham1 = alphas[0];
    alphas[0]  = rr_old / dAd;	// rr_old is taken from the last iteration respectively
    
   
    // r(k+1)
    assign_add_mul_r_32(r_up, Ad_up, (float) -alphas[0],N);
    assign_add_mul_r_32(r_dn, Ad_dn, (float) -alphas[0],N);

    // r(k+1)*r(k+1)
    rr_up  = square_norm_32(r_up, N, 1);
    rr_dn  = square_norm_32(r_dn, N, 1);
    rr     = rr_up + rr_dn;
    
      

    if((g_cart_id == 0) && (g_debug_level > 2)) printf("# CGMMSND_mixed: mms iteration j = %i: rr = %.6e\n", j, rr);

		 

    // aborting ?? // check wether precision is reached ...
    if ( ((check_abs)&&(rr <= eps_sq)) || ((check_rel)&&(rr <= eps_sq*r0r0)) ) 
    {
	if ((check_rel)&&(rr <= eps_sq*r0r0)) {
	  if((g_cart_id == 0) && (g_debug_level > 3)) printf("# CGMMSND_mixed: Reached relative solver precision of eps_rel = %.2e\n", eps_sq);
	}
      break;
   }
    
    // update alphas and zitas  
    // used later
    for(int im = 1; im < noshifts; im++) {
      gamma = zita[im]*alpham1/(alphas[0]*betas[0]*(1.-zita[im]/zitam1[im]) 
				+ alpham1*(1.+sigma[im]*alphas[0]));
      zitam1[im] = zita[im];
      zita[im] = gamma;
      alphas[im] = alphas[0]*zita[im]/zitam1[im];
    }  
    
    //check for reliable update
    res[0] = rr;
    for(int im=1; im<noshifts; im++) res[im] = rr * zita[im]; 
      
    rel_update = 0;
    for(int im = (noshifts-1); im >= 0; im--) {
      if( res[im] > maxres[im] ) maxres[im] = res[im];
      if( (res[im] < rel_delta*res0[im]) && (res0[im]<=maxres[im]) && (use_reliable) ) rel_update=1; 
      if( rel_update && ( trigger_shift == -1) ) trigger_shift = im;
    }     
    
    if(!rel_update)
    {
      // x_j(k+1) = x_j(k) + alpha_j*d_j(k) 
      // alphas are set above
      assign_add_mul_r_32(x_up, d_up, (float) alphas[0], N);   
      assign_add_mul_r_32(x_dn, d_dn, (float) alphas[0], N);
      
      
      for(int im = 1; im < noshifts; im++) {
	assign_add_mul_r_32(mms_x_up[im-1], mms_d_up[im-1], (float) alphas[im],  N);   
	assign_add_mul_r_32(mms_x_dn[im-1], mms_d_dn[im-1], (float) alphas[im],  N);  
      }  
   
      // beta = r(k+1)*r(k+1) / r(k)*r(k)
      betas[0] = rr / rr_old;
      rr_old = rr;  // for next iteration
      
      // d_0(k+1) = r(k+1) + beta*d_0(k) 
      assign_mul_add_r_32(d_up, (float) betas[0], r_up, N);  
      assign_mul_add_r_32(d_dn, (float) betas[0], r_dn, N); 
       
      // d_j(k+1) = zita*r(k+1) + beta*d_j(k)
      for(int im = 1; im < noshifts; im++) {
	betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
	assign_mul_add_mul_r_32(mms_d_up[im-1], r_up, (float) betas[im], (float) zita[im], N);
	assign_mul_add_mul_r_32(mms_d_dn[im-1], r_dn, (float) betas[im], (float) zita[im], N);
      }   
    }
    else{
      //reliable update
      if( (g_cart_id == 0) && (g_debug_level > 3) ){
	printf("# CGMMSND_mixed: Shift %d with offset squared %e triggered a reliable update\n", trigger_shift, sigma[trigger_shift]);
      }
      //add low prec solutions  
      assign_add_mul_r_32(x_up, d_up, (float) alphas[0], N); 
      assign_add_mul_r_32(x_dn, d_dn, (float) alphas[0], N); 
      
      addto_32(Pup[0], x_up, N);
      addto_32(Pdn[0], x_dn, N);	    
      for(int im = 1; im < noshifts; im++) {  
	assign_add_mul_r_32(mms_x_up[im-1], mms_d_up[im-1], alphas[im], N);
	assign_add_mul_r_32(mms_x_dn[im-1], mms_d_dn[im-1], alphas[im], N);	
	addto_32(Pup[im], mms_x_up[im-1], N);
        addto_32(Pdn[im], mms_x_dn[im-1], N);	
      }
      
      //add low precision for shift 0 only
      addto_32(x_up_d, x_up, N); 
      addto_32(x_dn_d, x_dn, N);      
 
      
      solver_params->M_ndpsi(Ax_up_d, Ax_dn_d, x_up_d,  x_dn_d);
      //add zero'th shift
      assign_add_mul_r(Ax_up_d, x_up_d, sigma[0], N);
      assign_add_mul_r(Ax_dn_d, x_dn_d, sigma[0], N);
      
      diff(r_up_d, Qup, Ax_up_d, N);         
      diff(r_dn_d, Qdn, Ax_dn_d, N); 
 
      rr_up = square_norm(r_up_d, N, 1);
      rr_dn = square_norm(r_dn_d, N, 1);
      rr    = rr_up + rr_dn;
      if ((g_cart_id == 0) && (g_debug_level > 3) ) printf("# CGMMSND_mixed: New residue after reliable update: %.6e\n", rr);
       
      //update res[im]
      res[0] = rr;

       
      if(res[trigger_shift] > res0[trigger_shift]){
	if(g_cart_id == 0) printf("# CGMMSND_mixed: Warning: residue of shift no %d got larger after rel. update\n", trigger_shift);
	//if this is the zero'th shift not getting better -> no further convergence, break
	if(trigger_shift == 0) break;
      }    
      
      //zero float fields
      zero_spinor_field_32(x_up, N);
      zero_spinor_field_32(x_dn, N);        
      for(int im = 1; im < noshifts; im++) {
	zero_spinor_field_32(mms_x_up[im-1], N);
	zero_spinor_field_32(mms_x_dn[im-1], N);  
      }
      
      //update the source
      assign_to_32(r_up, r_up_d, N);
      assign_to_32(r_dn, r_dn_d, N); 
      

      
      betas[0] = res[0]/rr_old;
      rr_old = rr;
      // d_0(k+1) = r(k+1) + beta*d_0(k)
      assign_mul_add_r_32(d_up, betas[0], r_up, N);
      assign_mul_add_r_32(d_dn, betas[0], r_dn, N);      
      // d_j(k+1) = r(k+1) + beta*d_j(k)
      for(int im = 1; im < noshifts; im++) {
	betas[im] = betas[0]*zita[im]*alphas[im]/(zitam1[im]*alphas[0]);
        assign_mul_add_mul_r_32(mms_d_up[im-1], r_up, (float) betas[im], (float) zita[im], N);
	assign_mul_add_mul_r_32(mms_d_dn[im-1], r_dn, (float) betas[im], (float) zita[im], N);
      } 
      
      //new maxres for the shift that initiated the reliable update
      res[trigger_shift] = res[0]*zita[trigger_shift]*zita[trigger_shift];
      res0[trigger_shift] = res[trigger_shift];  
      maxres[trigger_shift] = res[trigger_shift];
      trigger_shift = -1;
      no_rel_update ++;
    }	//reliable update	
    
    //check if some shift is converged
    for(int im = 1; im < noshifts; im++) {    
      if(j > 0 && (j % 10 == 0) && (im == noshifts-1)) {
	double sn = square_norm_32(mms_d_up[im-1], N, 1);
	sn +=       square_norm_32(mms_d_dn[im-1], N, 1);
	if(alphas[noshifts-1]*alphas[noshifts-1]*sn <= eps_sq) {
	  noshifts--;
	  if( (g_debug_level > 1) && (g_cart_id == 0) ) {
	    printf("# CGMMSND_mixed: at iteration %d removed one shift, %d remaining\n", j, noshifts);
	  }
	  //if removed we add the latest solution vector for this shift 	  
	  addto_32(Pup[im], mms_x_up[im-1], N);
          addto_32(Pdn[im], mms_x_dn[im-1], N);
	}
      }
    }
       
  }//LOOP
  
  if( (g_cart_id == 0) && (g_debug_level > 1) ) printf("Final mms residue: %.6e\n", rr);

  //add the latest solutions 
  for(int im = 0; im < noshifts; im++) {  
    if(im == 0){   
      addto_32(Pup[0], x_up, N);
      addto_32(Pdn[0], x_dn, N);        
    }
    else{     
      addto_32(Pup[im], mms_x_up[im-1], N);
      addto_32(Pdn[im], mms_x_dn[im-1], N);      
    }
  } 
  
  if(g_debug_level > 4){
    if(g_cart_id == 0) printf("# CGMMSND_mixed: Checking mms result:\n");
    //loop over all shifts (-> Nshift) 
    for(int im = 0; im < Nshift; im++){
      solver_params->M_ndpsi(sf[0], sf[1], Pup[im], Pdn[im]);
      assign_add_mul_r(sf[0], Pup[im] , shifts[im]*shifts[im], N);
      assign_add_mul_r(sf[1], Pdn[im] , shifts[im]*shifts[im], N);
      diff(sf[2], sf[0], Qup, N);
      diff(sf[3], sf[1], Qdn, N);
      rr_up = square_norm(sf[2], N, 1);
      rr_dn = square_norm(sf[3], N, 1);      
      rr = rr_up + rr_dn;
      if(g_cart_id == 0) printf("# CGMMSND_mixed: Shift[%d] squared residue: %e\n", im, rr);
    }
  }
  
 
  finalize_solver(sf, nr_sf);  
  finalize_solver_32(sf32, nr_sf32); 
 
  //free cg constants
  free(sigma); free(zitam1); free(zita); free(alphas); free(betas);    
  
  //free reliable update stuff
  free(res); free(res0); free(maxres);


  //if not converged -> return(-1)
  if(j<max_iter){
    return(j);
  }
  else{
    return(-1);
  }
}//










static unsigned int ini_mms_nd = 0;
static unsigned int nr_nd = 0;

static void init_mms_tm_nd_32(const unsigned int _nr, const unsigned int N) {
  if(ini_mms_nd == 0 || _nr > nr_nd) {
    if(nr_nd != 0) {
      free_mms_tm_nd_32();
    }
    nr_nd = _nr;

    x_up_qmms = (spinor32*)calloc(N*(nr_nd)+1,sizeof(spinor32));
    x_dn_qmms = (spinor32*)calloc(N*(nr_nd)+1,sizeof(spinor32));    
    d_up_qmms = (spinor32*)calloc(N*(nr_nd)+1,sizeof(spinor32));
    d_dn_qmms = (spinor32*)calloc(N*(nr_nd)+1,sizeof(spinor32));     
    mms_x_up = (spinor32**)calloc((nr_nd)+1,sizeof(spinor32*));
    mms_x_dn = (spinor32**)calloc((nr_nd)+1,sizeof(spinor32*));    
    mms_d_up = (spinor32**)calloc((nr_nd)+1,sizeof(spinor32*));
    mms_d_dn = (spinor32**)calloc((nr_nd)+1,sizeof(spinor32*));
    
    for(int i = 0; i < nr_nd; i++) {
      mms_x_up[i]=(spinor32*)(((unsigned long int)(x_up_qmms)+ALIGN_BASE32)&~ALIGN_BASE32) + i*N;
      mms_x_dn[i]=(spinor32*)(((unsigned long int)(x_dn_qmms)+ALIGN_BASE32)&~ALIGN_BASE32) + i*N;
      mms_d_up[i]=(spinor32*)(((unsigned long int)(d_up_qmms)+ALIGN_BASE32)&~ALIGN_BASE32) + i*N;
      mms_d_dn[i]=(spinor32*)(((unsigned long int)(d_dn_qmms)+ALIGN_BASE32)&~ALIGN_BASE32) + i*N;      
    }
    ini_mms_nd = 1;
  }
}

static void free_mms_tm_nd_32() {
  free(x_up_qmms); free(x_dn_qmms);
  free(d_up_qmms); free(d_dn_qmms);  
  free(mms_x_up); free(mms_x_dn);
  free(mms_d_up); free(mms_d_dn);  
  
  nr_nd = 0;
  ini_mms_nd = 0;
  return;
}




