/***********************************************************************
 *
 * Copyright (C) 2007 Craig McNeile
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
 *  Stout smear the gauge fields.
 *
 *  This routine overwrites the current copy of the gauge field.
 *
 *  The method is proposed in the paper:
 *   Analytic smearing of SU(3) link variables in lattice QCD.
 *   Colin Morningstar, Mike J. Peardon 
 *   Published in Phys.Rev.D69:054501,2004, hep-lat/0311018
 *
 *   Arguments:
 *             rho      : see equation 1 in paper above 
 *             no_iters : number of iterations of the algorithm
 *
 *
 * NOTES
 * -----
 * This version of stouting has been designed to get
 * the same answer as the stouting in chroma version
 * chroma-3.17.0.
 *
 * I have taken an empirical approach to getting various
 * signs correct. The occasional sign flips in the 
 * code were chosen to get agreement with the chroma 
 * code.
 *
 * written by Craig McNeile (1/2007)
 *
 **********************************************************************/

#include <lime.h>
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "get_staples.h"
#include "xchange_gauge.h"
#include "xchange.h"
#include "io/gauge.h"
#include "update_backward_gauge.h"
#include "stout_smear.h"

/*#include "stout_smear_aux_vars.h"*/

void  project_anti_herm(su3 *omega)  ;
void  print_su3(su3 *in)  ; 
void  check_su3_adj(su3 in, su3adj p) ;
void  scale_su3(su3 *in, double scale) ; 
su3   slow_expon(su3 in, int nterm)  ;
void  check_su3(su3 *in)  ;
void  print_config_to_screen(su3 **in);

int stout_smear_gauge_field(const double rho , const int no_iters) {
  const int dim=4 ; 
  int iter , mu , x; 
  su3 *gauge_wk[4] ; 
  su3  wk_staple  ; 
  su3 omega , Exp_p ; 
  su3adj p;
  su3  *gauge_local ; 
  su3  new_gauge_local ; 

  /*printf("Entering stout_smear_gauge_field\n");*/


  if(g_proc_id == 0 && g_debug_level > 3) {
    printf("DUMP OF g_gauge_field in STOUT\n");
    print_config_to_screen(g_gauge_field);

    printf("STOUT smearing the gauge fields\n") ; 
    printf("rho = %g number of iterations = %d\n",rho,no_iters) ; 
  }

  /* reserve memory  */
  for(mu = 0 ; mu < dim ; ++mu) {
    gauge_wk[mu] = calloc(VOLUME, sizeof(su3));
    if(errno == ENOMEM) {
      return(1);
    }
  }

  /* start of the the stout smearing **/

  for(iter = 0 ; iter < no_iters ; ++iter) {
    for(mu = 0 ; mu < dim  ; ++mu) {
      for(x= 0 ; x < VOLUME ; x++) {
	
	/*
	 *  we need to save all intermediate gauge configurations
	 *  because they are needed for the force back iteration in
	 *  "stout_smear_force.c"
	 */
	/*_su3_assign(g_gauge_field_smear_iterations[iter][x][mu], g_gauge_field[x][mu]);*/
	
	/* get staples */
	wk_staple = get_staples(x, mu, g_gauge_field) ; 
	scale_su3(&wk_staple, rho) ; 
	
	/* omega = staple * u^dagger */
	gauge_local = &g_gauge_field[x][mu];
	_su3_times_su3d(omega,wk_staple,*gauge_local);
	
	/* project out anti-hermitian traceless part */
	project_anti_herm(&omega) ; 
	
	/*  exponentiate */
	_trace_lambda(p,omega) ;
	/* -2.0 to get su3 to su3adjoint consistency ****/
	p.d1 /= -2.0 ; p.d2 /= -2.0 ; p.d3 /= -2.0 ; p.d4 /= -2.0 ; 
	p.d5 /= -2.0 ; p.d6 /= -2.0 ; p.d7 /= -2.0 ; p.d8 /= -2.0 ; 
	
	
	Exp_p = exposu3(p);
	
	/* new_gauge_local = Exp_p * gauge_local */
	_su3_times_su3(new_gauge_local,Exp_p,*gauge_local);
	gauge_wk[mu][x] = new_gauge_local ;  
	
      } /* end the loop over space-time */
    }
    /** update gauge field on this node **/
    for(mu = 0 ; mu < dim  ; ++mu) {
      for(x= 0 ; x < VOLUME ; ++x) {
	g_gauge_field[x][mu] = gauge_wk[mu][x] ;  
      }
    }
    if(g_debug_level > 3 && g_proc_id == 0) {
      printf("DUMP OF g_gauge_field in STOUT\n");
      print_config_to_screen(g_gauge_field);
    }
    
#ifdef MPI
    /** update boundaries for parallel stuff **/
    xchange_gauge();
#endif
    g_update_gauge_copy = 1;
    g_update_gauge_energy = 1;
    g_update_rectangle_energy = 1;
    
    /*
     *  here we save the intermediate smeares gauge fields a large array
     */
  } /* end loop over stout smearing iterations */
  
  /*    free up memory */
  for(mu=0 ; mu < dim ; ++mu) {
    free(gauge_wk[mu]);
  }
  
  if(g_debug_level > 3 && g_proc_id == 0) {
    printf("Leaving stout_smear_gauge_field\n");
  }
  return(0);
}

/*----------------------------------------------------------------------------*/

/*
   Project the hermitian part out of a
   su3 matrix. Also subtract the trace.
   Use the method in the Morninstar-Peardon
   paper, where I * anti-hermitian matrix
   is constructed.

   Note the sign flip at the end of the function.

 */

void  project_anti_herm(su3 *omega) {
  double tr_omega ; 

  tr_omega = (omega->c00.im + omega->c11.im + omega->c22.im) / 3.0 ; 

  omega->c00.re = 0.0 ;
  omega->c00.im = -omega->c00.im + tr_omega ;

  omega->c11.re = 0.0 ;
  omega->c11.im = -omega->c11.im + tr_omega ;

  omega->c22.re = 0.0 ;
  omega->c22.im = -omega->c22.im + tr_omega ;

  omega->c01.re -= omega->c10.re ;  omega->c01.re /= -2.0 ; 
  omega->c01.im += omega->c10.im ;  omega->c01.im /= -2.0 ; 
  omega->c10.re  = -omega->c01.re ; 
  omega->c10.im  =  omega->c01.im ; 

  omega->c02.re -= omega->c20.re ;  omega->c02.re /= -2.0 ; 
  omega->c02.im += omega->c20.im ;  omega->c02.im /= -2.0 ; 
  omega->c20.re  = -omega->c02.re ; 
  omega->c20.im  =  omega->c02.im ; 

  omega->c12.re -= omega->c21.re ;  omega->c12.re /= -2.0 ; 
  omega->c12.im += omega->c21.im ;  omega->c12.im /= -2.0 ; 
  omega->c21.re  = -omega->c12.re ; 
  omega->c21.im  =  omega->c12.im ; 


  /* flip sign  */
  scale_su3(omega, -1.0) ; 

}

/*##################### DEBUG #######################*/

void  load_config_from_file(su3 **in, char * filename) 
{
  int x, mu;
  su3 ** temp_su3;
  su3 * tsu3;
  temp_su3 = (su3**)calloc(VOLUME, sizeof(su3*));
  tsu3 = (su3*)calloc(4*VOLUME, sizeof(su3));
  temp_su3[0] = tsu3;
  for(x = 0; x < VOLUME; x++) {
    temp_su3[x] = temp_su3[x]+4;
  }

  for(x = 0; x < VOLUME; x++) {
    for(mu = 0; mu < 4; mu++) {
      _su3_assign(temp_su3[x][mu], g_gauge_field[x][mu]);
    }
  }
  read_gauge_field(filename);
  for(x = 0; x < VOLUME; x++) {
    for(mu = 0; mu < 4; mu++) {
      _su3_assign((in[x][mu]), g_gauge_field[x][mu]);
      _su3_assign(g_gauge_field[x][mu], temp_su3[x][mu]);
    }
  }
  free(temp_su3);
  free(tsu3);
}

/*----------------------------------------------------------------------------*/

void  print_config_to_screen(su3 **in) 
{
  int x, mu;
  for(x = 0; x < VOLUME; x++)
    for(mu = 0; mu < 4; mu++)
    {
      printf("x = %d  mu = %d\n", x, mu);
      /*print_su3_full_hex_precision(&(in[x][mu]));*/
      print_su3(&(in[x][mu]));
    }
}

/*----------------------------------------------------------------------------*/

void  print_scalar_complex_field_to_screen(complex *in_field) 
{
  int x;
  for(x = 0; x < VOLUME; x++)
    printf("x = %d:  %12.14f + I * %12.14f\n", x, in_field[x].re, in_field[x].im);
  printf("\n");
}

/*----------------------------------------------------------------------------*/

void  print_scalar_real_field_to_screen(double *in_field) 
{
  int x;
  for(x = 0; x < VOLUME; x++)
    printf("x = %d:  %12.14f\n", x, in_field[x]);
  printf("\n");
}

/*----------------------------------------------------------------------------*/

void  print_su3(su3 *in) {

  printf("[ %12.14f %12.14f , %12.14f %12.14f , %12.14f %12.14f  ] \n",
      in->c00.re,in->c00.im, 
      in->c01.re,in->c01.im, 
      in->c02.re,in->c02.im ) ; 
  printf("[ %12.14f %12.14f , %12.14f %12.14f , %12.14f %12.14f  ] \n",
      in->c10.re,in->c10.im, 
      in->c11.re,in->c11.im, 
      in->c12.re,in->c12.im ) ; 
  printf("[ %12.14f %12.14f , %12.14f %12.14f , %12.14f %12.14f  ] \n",
      in->c20.re,in->c20.im, 
      in->c21.re,in->c21.im, 
      in->c22.re,in->c22.im ) ; 


}

/*----------------------------------------------------------------------------*/

void  print_su3_full_hex_precision(su3 *in) {

  printf("[ %A %A , %A %A , %A %A  ] \n",
      in->c00.re,in->c00.im, 
      in->c01.re,in->c01.im, 
      in->c02.re,in->c02.im ) ; 
  printf("[ %A %A , %A %A , %A %A  ] \n",
      in->c10.re,in->c10.im, 
      in->c11.re,in->c11.im, 
      in->c12.re,in->c12.im ) ; 
  printf("[ %A %A , %A %A , %A %A  ] \n",
      in->c20.re,in->c20.im, 
      in->c21.re,in->c21.im, 
      in->c22.re,in->c22.im ) ; 


}

/*----------------------------------------------------------------------------*/

void  print_spinor(spinor *in) 
{
  printf("[ %f +i * %f,    %f + I * %f    ,    %f + I * %f  ] \n",
      in->s0.c0.re, in->s0.c0.im, in->s0.c1.re, in->s0.c1.im, in->s0.c2.re, in->s0.c2.im); 
  printf("[ %f +i * %f ,     %f + I * %f    ,    %f + I * %f  ] \n",
      in->s1.c0.re, in->s1.c0.im, in->s1.c1.re, in->s1.c1.im, in->s1.c2.re, in->s1.c2.im); 
  printf("[ %f +i * %f ,     %f + I * %f    ,    %f + I * %f  ] \n",
      in->s2.c0.re, in->s2.c0.im, in->s2.c1.re, in->s2.c1.im, in->s2.c2.re, in->s2.c2.im); 
  printf("[ %f +i * %f ,     %f + I * %f    ,    %f + I * %f  ] \n",
      in->s3.c0.re, in->s3.c0.im, in->s3.c1.re, in->s3.c1.im, in->s3.c2.re, in->s3.c2.im); 
}

void  print_su3_octave(su3 *in) 
{
  printf("[[(%12.14f)+i*(%12.14f) (%12.14f)+i*(%12.14f)  (%12.14f)+i*(%12.14f)]; [(%12.14f)+i*(%12.14f) (%12.14f)+i*(%12.14f)  (%12.14f)+i*(%12.14f)]; [(%12.14f)+i*(%12.14f) (%12.14f)+i*(%12.14f)  (%12.14f)+i*(%12.14f)]] \n",
      in->c00.re, in->c00.im, in->c01.re, in->c01.im, in->c02.re, in->c02.im, in->c10.re, in->c10.im, in->c11.re, in->c11.im, in->c12.re, in->c12.im, in->c20.re, in->c20.im, in->c21.re, in->c21.im, in->c22.re, in->c22.im ) ; 
}

/*##################### END OF DEBUG #######################*/

/*
   Check that a matrix is unitary
 */

void  check_su3(su3 *in) {

  su3 ans ; 
  _su3_times_su3d(ans,*in,*in); 

  printf("[ %f %f , %f %f , %f %f  ] \n",
      ans.c00.re,ans.c00.im, 
      ans.c01.re,ans.c01.im, 
      ans.c02.re,ans.c02.im ) ; 
  printf("[ %f %f , %f %f , %f %f  ] \n",
      ans.c10.re,ans.c10.im, 
      ans.c11.re,ans.c11.im, 
      ans.c12.re,ans.c12.im ) ; 
  printf("[ %f %f , %f %f , %f %f  ] \n",
      ans.c20.re,ans.c20.im, 
      ans.c21.re,ans.c21.im, 
      ans.c22.re,ans.c22.im ) ; 


}

/*----------------------------------------------------------------------------*/

/*
   Scale su3 matrix by a double.
   For performance this should be a MACRO
 */

void scale_su3(su3 *in, double scale) {

  in->c00.re *= scale  ; in->c00.im  *= scale  ;
  in->c01.re *= scale  ; in->c01.im  *= scale  ;
  in->c02.re *= scale  ; in->c02.im  *= scale  ;

  in->c10.re *= scale  ; in->c10.im *= scale  ;
  in->c11.re *= scale  ; in->c11.im *= scale  ;
  in->c12.re *= scale  ; in->c12.im *= scale  ;

  in->c20.re *= scale  ; in->c20.im  *= scale  ;
  in->c21.re *= scale  ; in->c21.im *= scale  ;
  in->c22.re *= scale  ; in->c22.im  *= scale  ;


}

/*----------------------------------------------------------------------------*/

/*
   Check the conversion of su3 adjoint 
   to su3
 */

void check_su3_adj(su3 in, su3adj p) {
  su3 out ;

  printf("INPUT matrix \n") ; 
  print_su3(&in); 

  printf("lambda_[1 - 4] =  %f %f %f %f \n",p.d1,p.d2, p.d3,p.d4) ; 
  printf("lambda_[5 - 8] =  %f %f %f %f \n",p.d5,p.d6, p.d7,p.d8) ; 

  /**
    p.d1 /= 2.0 ; p.d2 /= 2.0 ; p.d3 /= 2.0 ; p.d4 /= 2.0 ; 
    p.d5 /= 2.0 ; p.d6 /= 2.0 ; p.d7 /= 2.0 ; p.d8 /= 2.0 ; 
   **/

  _make_su3(out,p) ;

  printf("OUTPUT matrix \n") ; 
  print_su3(&out); 

}

/*----------------------------------------------------------------------------*/

/*
   Power series exponentiation.
   This does not involve the su3adj stuff.

   There is another power series routine
   in the expo.c file: exposu3_check(p, 20) ; 

 */

su3 slow_expon(su3 in, int nterm)  {
  su3 out , out_new , aa , aa_old , aa_scale ;
  int i ; 
  double fact = 1 ; 

  _su3_one(out) ; 
  aa = in ; 
  aa_old = aa ; 

  _su3_plus_su3(out_new,out,aa) ; 
  out = out_new ;  

  for(i=2 ; i < nterm ; ++i) {

    _su3_times_su3(aa,in,aa_old);
    aa_old = aa ; 
    aa_scale = aa ; 
    fact *= i ; 
    scale_su3(&aa_scale, 1.0/fact ) ; 

    _su3_plus_su3(out_new,out,aa_scale) ; 
    out = out_new ;  
  }

  return(out); 
}

