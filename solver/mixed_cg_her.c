/***********************************************************************
 * Copyright (C) 2013 Florian Burger
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
 * 
 *  
 * File: mixed_cg_her.c
 *
 * CG solver for hermitian f only!
 *
 * The externally accessible functions are
 *
 *
 *   int cg(spinor * const P, spinor * const Q, double m, const int subtract_ev)
 *     CG solver
 *
 * input:
 *   m: Mass to be use in D_psi
 *   subtrac_ev: if set to 1, the lowest eigenvectors of Q^2 will
 *               be projected out.
 *   Q: source
 * inout:
 *   P: initial guess and result
 * 
 *
 **************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "su3.h"
#include "linalg_eo.h"
#include "start.h"
#include "operator/tm_operators_32.h"
#include "solver/matrix_mult_typedef.h"

#include "solver_field.h"
#include "solver/mixed_cg_her.h"
#include "gettime.h"


#include "bgq.h"
#include "bgq2.h"
#include "operator/halfspinor_hopping32.h"
void test(){
 
 su3_32 ALIGN32 Utest[2];

 Utest[0].c00 =  0.0f+0.0fI;
 Utest[0].c01 =  0.0f+0.0fI;
 Utest[0].c02 =  0.0f+0.0fI; 
 Utest[0].c10 =  0.0f+0.0fI;
 Utest[0].c11 =  0.0f+0.0fI;
 Utest[0].c12 =  0.0f+0.0fI; 
 Utest[0].c20 =  0.0f+0.0fI;
 Utest[0].c21 =  0.0f+0.0fI;
 Utest[0].c22 =  0.0f+0.0fI; 



 Utest[1].c00 =  2.0f+0.0fI;
 Utest[1].c01 =  0.0f+0.0fI;
 Utest[1].c02 =  0.0f+0.0fI; 
 Utest[1].c10 =  0.0f+0.0fI;
 Utest[1].c11 =  2.0f+0.0fI;
 Utest[1].c12 =  0.0f+0.0fI; 
 Utest[1].c20 =  0.0f+0.0fI;
 Utest[1].c21 =  0.0f+0.0fI;
 Utest[1].c22 =  2.0f+0.0fI; 


  
  spinor32 ALIGN32 a;
  a.s0.c0 = 1.0f+0.0fI; 
  a.s0.c1 = 1.0f+0.0fI; 
  a.s0.c2 = 1.0f+0.0fI; 
   
  a.s1.c0 = 2.0f+0.0fI; 
  a.s1.c1 = 2.0f+0.0fI; 
  a.s1.c2 = 2.0f+0.0fI; 
  
  a.s2.c0 = 3.0f+0.0fI; 
  a.s2.c1 = 3.0f+0.0fI; 
  a.s2.c2 = 3.0f+0.0fI; 
  
  a.s3.c0 = 4.0f+0.0fI; 
  a.s3.c1 = 4.0f+0.0fI; 
  a.s3.c2 = 4.0f+0.0fI; 
   

  spinor32 ALIGN32 out;
  out.s0.c0 = 11.0f;
  out.s0.c1 = 11.0f;
  out.s0.c2 = 11.0f;
  out.s1.c0 = 11.0f;
  out.s1.c1 = 11.0f;
  out.s1.c2 = 11.0f;

  _declare_hregs();
  _Complex float ALIGN32 ka0_32 = 1.0f;
  spinor32 * s ALIGN32;
  su3_32 * U ALIGN32;
  halfspinor32 * phi2 ALIGN32;
  
  U=&(Utest[1]);
  s=&a;
  phi2=(halfspinor32*)&(out);

#if ( defined BGQ && defined XLC)
 #define _prefetch_spinor_32
 #define _prefetch_su3_32 
  printf("Testing qpx...\n");
  /*
  _vec_load2_32(rs0, rs1, rs2, s->s0);				
  _vec_load2_32(rs3, rs4, rs5, s->s1);				
  _vec_load2_32(rs6, rs7, rs8, s->s2);				
  _vec_load2_32(rs9, rs10, rs11, s->s3);			
  _vec_add_to2(r0, r1, r2, rs0, rs1, rs2, rs6, rs7, rs8);	
  _vec_add_to2(r3, r4, r5, rs3, rs4, rs5, rs9, rs10, rs11);
  rtmp = vec_ld2a(0L, (float*) &ka0_32);			
  _vec_su3_multiply_double2_32(U);		
  _vec_cmplx_mul_double2(r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, rtmp); 
  _vec_store2_32(out.s0, r0, r1, r2);
  _vec_store2_32(out.s1, r3, r4, r5);
  */
  _vec_load_32(rs0, rs1, s->s0);		
  _vec_load16_32(rs2, rs3, s->s1, rtmp);
  _vec_load_32(rs4, rs5, s->s2);
  _vec_load16_32(rs6, rs7, s->s3, rtmp);
  _vec_add(r0, r1, rs0, rs1, rs4, rs5);	
  _vec_add(r2, r3, rs2, rs3, rs6, rs7);
  _vec_su3_multiply_double2c_32(U);
  rtmp = vec_ld2(0, (float*) &ka0_32);
  _vec_cmplx_mul_double2c(r0, r1, r2, r4, r5, r6, rtmp);
  _vec_store_halfspinor_32(phi2->s0, r0, r1, r2);

#endif

  printf("IN:  (%f,%f) (%f,%f) (%f,%f)\n", creal(a.s0.c0), cimag(a.s0.c0), creal(a.s0.c1), cimag(a.s0.c1), creal(a.s0.c2), cimag(a.s0.c2));
  printf("OUT: (%f,%f) (%f,%f) (%f,%f)\n", creal(out.s0.c0), cimag(out.s0.c0), creal(out.s0.c1), cimag(out.s0.c1), creal(out.s0.c2), cimag(out.s0.c2));

//exit(10);

return;
}



/* P output = solution , Q input = source */
int mixed_cg_her(spinor * const P, spinor * const Q, const int max_iter, 
		 double eps_sq, const int rel_prec, const int N, matrix_mult f, matrix_mult32 f32) {

  int i = 0, iter = 0, j = 0;
  float sqnrm = 0., sqnrm2, squarenorm;
  float pro, err, alpha_cg, beta_cg;
  double sourcesquarenorm, sqnrm_d, squarenorm_d;
  spinor *delta, *y, *xhigh;
  spinor32 *x, *stmp;
  spinor ** solver_field = NULL;
  spinor32 ** solver_field32 = NULL;  
  const int nr_sf = 3;
  const int nr_sf32 = 4;

  int max_inner_it = 2800;
  int N_outer = max_iter/max_inner_it;

  //test();

  if(N == VOLUME) {
    init_solver_field(&solver_field, VOLUMEPLUSRAND, nr_sf);    
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND, nr_sf32);
  }
  else {
    init_solver_field(&solver_field, VOLUMEPLUSRAND/2, nr_sf);
    init_solver_field32(&solver_field32, VOLUMEPLUSRAND/2, nr_sf32);    
  }

/*
misaligned = ((long) Q) & 0x1F ;
if(misaligned){
    printf("x is misaligned!!\n");
}
  else{
    printf("x correctly aligned!!\n");
}
*/
  squarenorm_d = square_norm(Q, N, 1);
  //printf("sqarenorm_d = %f\n", squarenorm_d);
  sourcesquarenorm = squarenorm_d;
  sqnrm_d = squarenorm_d;
  
  //maybe that is wrong
  squarenorm_d = square_norm(P, N, 1);
 
  //printf("sqarenorm_d = %f\n", squarenorm_d);
 
  delta = solver_field[0];
  y = solver_field[1];
  xhigh = solver_field[2];
  x = solver_field32[3];  
  
  assign(delta, Q, N);


//printf("In mixed_cg solver...\n");
//printf("volume is: %d\n",N); 
 /* here comes a small test
   spinor32 * help_low = solver_field32[0];
   spinor * help_high = solver_field[0];
   assign_to_32(help_low, Q, N);
   assign(help_high, Q, N);
   printf("square_norm(Q_high) = %e\n", square_norm(help_high,N,1));
   printf("square_norm(Q_low) = %e\n", square_norm_32(help_low,N,1));  
   f32(solver_field32[1], help_low);
   f(solver_field[1], help_high);
   
   assign_to_64(solver_field[2], solver_field32[1], N);
   diff(solver_field[3], solver_field[1], solver_field[2], N);
   sqnrm = square_norm(solver_field[3], N, 1);
   printf("Operator 32 test: (square_norm) / (spinor component) = %.8e\n", sqnrm/24.0/VOLUME);
  // exit(1);
  /* end of small test*/
  
  
  /* small performance test */
  
//   int Nhit = 40;
//   double t1,t2,dt,sdt;
//   double antioptaway=0.0;
//   t1 = gettime();
//   antioptaway=0.0;
//   for (i=0;i<Nhit;i++) {
//      f(y, delta);
//      f(delta,y);
//   }
//      antioptaway+=creal(delta[0].s0.c0);  
//   t2 = gettime();
//   dt = t2-t1;
//   sdt=2*1.0e6f*dt/((double)(Nhit*(N)));
//   printf("antioptaway = %e\n", antioptaway);
//   printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt), (int)sizeof(spinor)/3);
//   assign(delta, Q, N);
//   
//   
//   assign_to_32(x, delta, N);
//   float antioptaway_f=0.0;
//   t1 = gettime();
//   antioptaway_f=0.0;
//   for (i=0;i<Nhit;i++) {
//      f32(solver_field32[0], x);
//      f32(x,solver_field32[0]);
// 
//   }
//   antioptaway_f+=creal(x[0].s0.c0);
//   t2 = gettime();
//   dt = t2-t1;
//   sdt=2*1.0e6f*dt/((double)(Nhit*(N)));
//   printf("antioptaway = %e\n", antioptaway_f);
//   printf("# Communication switched on:\n# (%d Mflops [%d bit arithmetic])\n", (int)(1608.0f/sdt), (int)sizeof(spinor32)/3);
//   
  /* end of small performance test */
  
  
  if(squarenorm_d > 1.e-7) { 
    /* if a starting solution vector different from zero is chosen */ 
    printf("We have a non-zero starting solution -> using it\n");
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    printf("sqnrm_d = %f\n", sqnrm_d);
    if(((sqnrm_d <= eps_sq) && (rel_prec == 0)) || ((sqnrm_d <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      finalize_solver32(solver_field32, nr_sf32);      
      return(0);
    }
  }
  

  for(i = 0; i < N_outer; i++) {

    //g_sloppy_precision = 1;
    //g_sloppy_precision_flag = 1;
    /* main CG loop in lower precision */
    zero_spinor_field_32(x, N);
    zero_spinor_field_32(solver_field32[0], N);   
    assign_to_32(solver_field32[1], delta, N);
    assign_to_32(solver_field32[2], delta, N);
    
    sqnrm = (float) sqnrm_d;
    sqnrm2 = sqnrm;
    for(j = 0; j <= max_inner_it; j++) {
      f32(solver_field32[0], solver_field32[2]); 
      pro = scalar_prod_r_32(solver_field32[2], solver_field32[0], N, 1);
      alpha_cg = sqnrm2 / pro;
      
      assign_add_mul_r_32(x, solver_field32[2], alpha_cg, N);
      
      assign_mul_add_r_32(solver_field32[0], -alpha_cg, solver_field32[1], N);      
      
      err = square_norm_32(solver_field32[0], N, 1);

      if(g_proc_id == g_stdio_proc && g_debug_level > 2) {
	printf("inner CG: %d res^2 %g\n", iter+j, err);
	fflush(stdout);
      }
    
      //if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
      if((err <= 6.0e-6*sqnrm)|| (j==max_inner_it) ||  ((1.3*err <= eps_sq) && (rel_prec == 0)) || ((1.3*err <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
	break;
      }
      beta_cg = err / sqnrm2;
      assign_mul_add_r_32(solver_field32[2], beta_cg, solver_field32[0], N);
      stmp = solver_field32[0];
      solver_field32[0] = solver_field32[1];
      solver_field32[1] = stmp;
      sqnrm2 = err;
      
    }
    /* end main CG loop */
    iter += j;
    //g_sloppy_precision = 0;
    //g_sloppy_precision_flag = 0;
    assign_to_64(xhigh, x, N);    
    
//     f(y, xhigh);
//     diff(P, delta, y, N);
//     sqnrm = square_norm(P, N, 1);
//     printf("mixed CG: true residue %d\t%g\t\n",iter, sqnrm);   
//     exit(1);
    
    
    add(P, P, xhigh, N);
    f(y, P);
    diff(delta, Q, y, N);
    sqnrm_d = square_norm(delta, N, 1);
    if(g_debug_level > 0 && g_proc_id == g_stdio_proc) {
      printf("mixed CG: last inner residue: %g\t\n", err);
      printf("mixed CG: true residue %d %g\t\n",iter, sqnrm_d); fflush(stdout);
    }

    if(((sqnrm_d <= eps_sq) && (rel_prec == 0)) || ((sqnrm_d <= eps_sq*sourcesquarenorm) && (rel_prec == 1))) {
      finalize_solver(solver_field, nr_sf);
      finalize_solver32(solver_field32, nr_sf32); 
      return(iter+i);
    }
    iter++;
  }
  finalize_solver(solver_field, nr_sf);
  finalize_solver32(solver_field32, nr_sf32); 
  return(-1);
}







