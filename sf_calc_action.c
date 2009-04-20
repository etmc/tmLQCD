/*******************************************
*
* FILE: calc_action.c
*
* Author: Jenifer Gonzalez Lopez
*
********************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "geometry_eo.h"
#include "sf_calc_action.h"

/**************************************************************************************************/

/* the next function imposes Dirichlet b.c.
by just setting all the gauge links in the time direction (from t on) to zero.
Note that the rest of the links at the boundaries (spatial links) are not yet touched here */
void dirichlet_boundary_conditions(int t) {
  
  int ix;
  
  for (ix=0;ix<VOLUME;ix++){
 
    if (g_t[ix] == t) {

      _su3_zero(g_gauge_field[ix][0]);

    }   
  }
}

/* this function sets all the spatial links located at the time boundaries to one
   it does nothing to the time-like links */
void dirichlet_boundary_conditions_spatial_links_to_one(int t) {
  
  int ix;
  
  for (ix=0;ix<VOLUME;ix++){
 
    if (g_t[ix] == 0 || g_t[ix] == t) {

      _su3_one(g_gauge_field[ix][1]);
      _su3_one(g_gauge_field[ix][2]);
      _su3_one(g_gauge_field[ix][3]);

    }   
  }
}


/* this function sets all the links in the lattice to one */
void set_all_links_to_one() {
  
  int ix;
  
  for (ix=0;ix<VOLUME;ix++){
 
      _su3_one(g_gauge_field[ix][0]);
      _su3_one(g_gauge_field[ix][1]);
      _su3_one(g_gauge_field[ix][2]);
      _su3_one(g_gauge_field[ix][3]);

  }
}

/* this function sets all the links in the lattice to one, except the temporal links at t wich are set to zero */
void set_all_links_to_one_with_dirichlet(int t) {
  
  int ix;
  
  for (ix=0;ix<VOLUME;ix++){
    
    if (g_t[ix] == t) {
      
      _su3_zero(g_gauge_field[ix][0]);
      _su3_one(g_gauge_field[ix][1]);
      _su3_one(g_gauge_field[ix][2]);
      _su3_one(g_gauge_field[ix][3]);
      
    } 
    
    else {
      
      _su3_one(g_gauge_field[ix][0]);
      _su3_one(g_gauge_field[ix][1]);
      _su3_one(g_gauge_field[ix][2]);
      _su3_one(g_gauge_field[ix][3]);  
    }
    
    
  }
}


/* it calculates an su3 matrix "u" which is gonna be the (lattice) spatially constant abelian field */
#define _su3_spatially_constant_abelian_field(u, p1, p2, p3)	\
  (u).c00.re = cos(p1);						\
  (u).c00.im = sin(p1);						\
  (u).c01.re = 0.0;						\
  (u).c01.im = 0.0;						\
  (u).c02.re = 0.0;						\
  (u).c02.im = 0.0;						\
  (u).c10.re = 0.0;						\
  (u).c10.im = 0.0;						\
  (u).c11.re = cos(p2);						\
  (u).c11.im = sin(p2);						\
  (u).c12.re = 0.0;						\
  (u).c12.im = 0.0;						\
  (u).c20.re = 0.0;						\
  (u).c20.im = 0.0;						\
  (u).c21.re = 0.0;						\
  (u).c21.im = 0.0;						\
  (u).c22.re = cos(p3);						\
  (u).c22.im = sin(p3);


/* it calculates an su3 matrix "u" which is gonna be the (continuum) spatially constant abelian field */
#define _su3_spatially_constant_abelian_field_continuum(u, p1, p2, p3)	\
  (u).c00.re = 0.0;							\
  (u).c00.im = p1;							\
  (u).c01.re = 0.0;							\
  (u).c01.im = 0.0;							\
  (u).c02.re = 0.0;							\
  (u).c02.im = 0.0;							\
  (u).c10.re = 0.0;							\
  (u).c10.im = 0.0;							\
  (u).c11.re = 0.0;							\
  (u).c11.im = p2;							\
  (u).c12.re = 0.0;							\
  (u).c12.im = 0.0;							\
  (u).c20.re = 0.0;							\
  (u).c20.im = 0.0;							\
  (u).c21.re = 0.0;							\
  (u).c21.im = 0.0;							\
  (u).c22.re = 0.0;							\
  (u).c22.im = p3;


/* it just prints on the screen the su3 matrix "u" */
void print_su3_matrix (su3 u) {
  
  printf(" %e i %e  %e i %e  %e i %e \n",
	 (u).c00.re,  (u).c00.im,  (u).c01.re,  (u).c01.im,  (u).c02.re,  (u).c02.im);
  printf(" %e i %e  %e i %e  %e i %e \n",
	 (u).c10.re,  (u).c10.im,  (u).c11.re,  (u).c11.im,  (u).c12.re,  (u).c12.im);
  printf(" %e i %e  %e i %e  %e i %e \n",
	 (u).c20.re,  (u).c20.im,  (u).c21.re,  (u).c21.im,  (u).c22.re,  (u).c22.im);
  printf("\n");
}

/* this function gives us the boundary gauge field: spatially constant abelian field */
void sf_boundary_conditions_spatially_constant_abelian_field(int t, double eta) {
  
  int ix;
  double pi;
  double phi1_0, phi2_0, phi3_0;
  double phi1_T, phi2_T, phi3_T;

  pi = acos(-1.);

  phi1_0 = eta - pi/3.0;
  phi2_0 = - 0.5 * eta;
  phi3_0 = - 0.5 * eta + pi/3.0;
  
  phi1_T = - phi1_0 - (4.0*pi)/3.0;
  phi2_T = - phi3_0 + (2.0*pi)/3.0;
  phi3_T = - phi2_0 + (2.0*pi)/3.0;

  phi1_0 /= (double)LX; 
  phi2_0 /= (double)LX; 
  phi3_0 /= (double)LX; 
  
  phi1_T /= (double)LX; 
  phi2_T /= (double)LX; 
  phi3_T /= (double)LX; 

  for (ix=0;ix<VOLUME;ix++){
    
    if (g_t[ix] == 0) {
      
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][1], phi1_0, phi2_0, phi3_0);
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][2], phi1_0, phi2_0, phi3_0);
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][3], phi1_0, phi2_0, phi3_0);

    } 

    
    if (g_t[ix] == t) {
  
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][1], phi1_T, phi2_T, phi3_T);
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][2], phi1_T, phi2_T, phi3_T);
      _su3_spatially_constant_abelian_field(g_gauge_field[ix][3], phi1_T, phi2_T, phi3_T);

    }   


  }
  
}



/*** MEASUREMENTS ***/

/* it calculates the plaquette (see notes for notation) for PBC */
double measure_plaquette() {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<3;mu1++){
      
      ix1=g_iup[ix][mu1];
      
      for (mu2=mu1+1;mu2<4;mu2++){ 
	
	ix2=g_iup[ix][mu2];
	
	v=&g_gauge_field[ix][mu1];
	w=&g_gauge_field[ix1][mu2];
	
	_su3_times_su3(pr1,*v,*w);
	
	
	v=&g_gauge_field[ix][mu2];
	w=&g_gauge_field[ix2][mu1];
	
	_su3_times_su3(pr2,*v,*w);
	
	
	_trace_su3_times_su3d(ac,pr1,pr2);
	
	sum += ac;
	
      }
    }
    
  }
  
  ga = sum*2.0;

  return ga;
  
}


/* it calculates the rectangle (see notes for notation) for PBC */
double measure_rectangle() {

  int ix,ix1,ix2,mu1,mu2;
  int ix12,ix22;                           /* My changes */
  static su3 pr_r1,pr_r2,pr_r11,pr_r22;    /* My changes */
  static double ga,ac; 
  su3 *v_r1,*w_r1,*z_r1;                   /* My changes */
  su3 *v_r2,*w_r2,*z_r2;                   /* My changes */
  double sum = 0;                          /* My changes */


  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<4;mu1++){
      
      ix1 = g_iup[ix][mu1];
      
      for (mu2=0;mu2<4;mu2++){
	if( mu1 != mu2 ){
	  
	  ix2  = g_iup[ix][mu2];
	  ix12 = g_iup[ix1][mu2];
	  ix22 = g_iup[ix2][mu2];
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu2];
	  z_r1 = &g_gauge_field[ix12][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu2];
	  z_r2 = &g_gauge_field[ix22][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr_r22,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr_r11,pr_r22);
	  
	  sum += ac;
	}
      }
    }
  }
  
  ga = sum*2.0;

  return ga;

}

/* SF functions: */

/* it calculates the plaquette for SF b.c. without any improvement coefficient */
/* "hard-coded": boundaries and bulk defined in the same function */
double measure_plaquette_sf_weights(int t) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<3;mu1++){
      
      ix1=g_iup[ix][mu1];
      
      for (mu2=mu1+1;mu2<4;mu2++){ 
	
	ix2=g_iup[ix][mu2];
	
	v=&g_gauge_field[ix][mu1];
	w=&g_gauge_field[ix1][mu2];
	
	_su3_times_su3(pr1,*v,*w);
	
	
	v=&g_gauge_field[ix][mu2];
	w=&g_gauge_field[ix2][mu1];
	
	_su3_times_su3(pr2,*v,*w);
	
	
	_trace_su3_times_su3d(ac,pr1,pr2);

	if (g_t[ix] == 0) {
	  
	  if (mu1 != 0 && mu2 != 0) {
	    
	    ac *= 0.5;
	    
	  }
	  
	} 
	
	
	if (g_t[ix] == t) {
	  
	  if (mu1 != 0 && mu2 != 0) {
	    
	    ac *= 0.5;
	    
	  }	  
	  
	} 
	
	sum += ac;
	
      }
    }
    
  }
  
  ga = sum*2.0;

  return ga;
  
}


/* it calculates the plaquette for SF b.c. WITH the (alpha collab.) improvement coefficients */
/* "hard-coded": boundaries and bulk defined in the same function */
double measure_plaquette_sf_weights_improvement(int t, double cs, double ct) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;


  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<3;mu1++){
      
      ix1=g_iup[ix][mu1];
      
      for (mu2=mu1+1;mu2<4;mu2++){ 
	
	ix2=g_iup[ix][mu2];
	
	v=&g_gauge_field[ix][mu1];
	w=&g_gauge_field[ix1][mu2];
	
	_su3_times_su3(pr1,*v,*w);
	
	
	v=&g_gauge_field[ix][mu2];
	w=&g_gauge_field[ix2][mu1];
	
	_su3_times_su3(pr2,*v,*w);
	
	
	_trace_su3_times_su3d(ac,pr1,pr2);

	if (g_t[ix] == 0) {
	  
	  if (mu1 != 0 && mu2 != 0) {

	    ac *= cs;

	  }

	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {

	    ac *= ct;

	  }
	
	} 
	
	
	if (g_t[ix] == t) {

	  if (mu1 != 0 && mu2 != 0) {
	    
	    ac *= cs;
	    
	  }

	}

	if (g_t[ix] == (t-1)) {
	  
	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	    
	    ac *= ct;
	    
	  }	  
	  
	} 
	
	sum += ac;
	
      }
    }
    
  }
  
  ga = sum*2.0;

  return ga;
  
}

/*-----------------------------------------------------------------------------------------------------------*/

/* the next three functions summed up will give:
the plaquette for SF b.c. withOUT the improvement coefficients */
/* differently as in the previous cases, now bulk and boundaries are calculated in different functions */

/* (1) bulk: */
double measure_plaquette_sf_weights_bulk(int t) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  for (x0 = 1 ; x0 < t ; x0++) {
    
    for (x1 = 0 ; x1 < LX ; x1++) {
      
      for (x2 = 0 ; x2 < LY ; x2++) {
	
	for (x3 = 0 ; x3 < LZ ; x3++) {
	  
	  ix = Index(x0, x1, x2, x3);
	  
	  for (mu1=0;mu1<3;mu1++) {
	    
	    ix1=g_iup[ix][mu1];
	    
	    for (mu2=mu1+1;mu2<4;mu2++) { 
	      
	      ix2=g_iup[ix][mu2];
	      
	      v=&g_gauge_field[ix][mu1];
	      w=&g_gauge_field[ix1][mu2];
	      
	      _su3_times_su3(pr1,*v,*w);
	      
	      
	      v=&g_gauge_field[ix][mu2];
	      w=&g_gauge_field[ix2][mu1];
	      
	      _su3_times_su3(pr2,*v,*w);
	      
	      
	      _trace_su3_times_su3d(ac,pr1,pr2);
	      
	      sum += ac;
	      
	    }
	  }	      
	  
	}
      }
    }
  }
  
  ga = sum*2.0;
  
  return ga;
  
}

/* (2) boundary at 0 */
double measure_plaquette_sf_weights_boundary_0 () {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  x0 = 0;
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=0;mu1<3;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  
	  for (mu2=mu1+1;mu2<4;mu2++) { 
	    
	    ix2=g_iup[ix][mu2];
	    
	    v=&g_gauge_field[ix][mu1];
	    w=&g_gauge_field[ix1][mu2];
	    
	    _su3_times_su3(pr1,*v,*w);
	    
	    
	    v=&g_gauge_field[ix][mu2];
	    w=&g_gauge_field[ix2][mu1];
	    
	    _su3_times_su3(pr2,*v,*w);
	    
	    
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    if (mu1 != 0 && mu2 != 0) {
	      
	      ac *= 0.5;
	      
	    }
	    
	    sum += ac;
	    
	  }
	}	      	
	
      }
    }
  }

  ga = sum*2.0;
  
  return ga;
  
}

/* (3) boundary at t */
double measure_plaquette_sf_weights_boundary_t (int t) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  x0 = t;
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=0;mu1<3;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  
	  for (mu2=mu1+1;mu2<4;mu2++) { 
	    
	    ix2=g_iup[ix][mu2];
	    
	    v=&g_gauge_field[ix][mu1];
	    w=&g_gauge_field[ix1][mu2];
	    
	    _su3_times_su3(pr1,*v,*w);
	    
	    
	    v=&g_gauge_field[ix][mu2];
	    w=&g_gauge_field[ix2][mu1];
	    
	    _su3_times_su3(pr2,*v,*w);
	    
	    
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    if (mu1 != 0 && mu2 != 0) {
	      
	      ac *= 0.5;
	      
	    }
	    
	    sum += ac;
	    
	  }
	}	      	
	
      }
    }
  }
  
  ga = sum*2.0;
  
  return ga;
  
}


/*-------------------------------------------------------------------------------------------------------------*/



/* the next four functions summed up will give:
the plaquette for SF b.c. WITH the improvement coefficients (alpha collab.) */
/* differently as in the previous cases, now bulk and boundaries are calculated in different functions */

/* (1) bulk: */
double measure_plaquette_sf_weights_improved_bulk(int t) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  for (x0 = 1 ; x0 < (t-1) ; x0++) {
    
    for (x1 = 0 ; x1 < LX ; x1++) {
      
      for (x2 = 0 ; x2 < LY ; x2++) {
	
	for (x3 = 0 ; x3 < LZ ; x3++) {
	  
	  ix = Index(x0, x1, x2, x3);
	  
	  for (mu1=0;mu1<3;mu1++) {
	    
	    ix1=g_iup[ix][mu1];
	    
	    for (mu2=mu1+1;mu2<4;mu2++) { 
	      
	      ix2=g_iup[ix][mu2];
	      
	      v=&g_gauge_field[ix][mu1];
	      w=&g_gauge_field[ix1][mu2];
	      
	      _su3_times_su3(pr1,*v,*w);
	      
	      
	      v=&g_gauge_field[ix][mu2];
	      w=&g_gauge_field[ix2][mu1];
	      
	      _su3_times_su3(pr2,*v,*w);
	      
	      
	      _trace_su3_times_su3d(ac,pr1,pr2);
	      
	      sum += ac;
	      
	    }
	  }	      
	  
	}
      }
    }
  }
  
  ga = sum*2.0;
  
  return ga;
  
}

/* (2) boundary at 0 */
double measure_plaquette_sf_weights_improved_boundary_0 (double cs, double ct) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  x0 = 0;


  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);


	for (mu1=0;mu1<3;mu1++){
	  
	  ix1=g_iup[ix][mu1];
	  
	  for (mu2=mu1+1;mu2<4;mu2++){ 
	    
	    ix2=g_iup[ix][mu2];
	    
	    v=&g_gauge_field[ix][mu1];
	    w=&g_gauge_field[ix1][mu2];
	    
	    _su3_times_su3(pr1,*v,*w);
	    
	    
	    v=&g_gauge_field[ix][mu2];
	    w=&g_gauge_field[ix2][mu1];
	    
	    _su3_times_su3(pr2,*v,*w);
	    
	
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    
	    if (mu1 != 0 && mu2 != 0) {
	      
	      ac *= cs;
	      
	    }
	    
	    if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	      
	      ac *= ct;
	      
	    }
	    
	    sum += ac;
	    
	  }
	}
	
      }
    }
  }
  
  ga = sum*2.0;

  return ga;  

}

/*(3) boundary at t */
double measure_plaquette_sf_weights_improved_boundary_t (int t, double cs) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  x0 = t;
  
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	
	for (mu1=0;mu1<3;mu1++){
	  
	  ix1=g_iup[ix][mu1];
	  
	  for (mu2=mu1+1;mu2<4;mu2++){ 
	    
	    ix2=g_iup[ix][mu2];
	    
	    v=&g_gauge_field[ix][mu1];
	    w=&g_gauge_field[ix1][mu2];
	    
	    _su3_times_su3(pr1,*v,*w);
	    
	    
	    v=&g_gauge_field[ix][mu2];
	    w=&g_gauge_field[ix2][mu1];
	    
	    _su3_times_su3(pr2,*v,*w);
	    
	    
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    
	    if (mu1 != 0 && mu2 != 0) {
	      
	      ac *= cs;
	      
	    }
	    
	    sum += ac;
	    
	  }
	}
	
      }
    }
  }
  
  ga = sum*2.0;

  return ga;
  
}

/* (4) boundary at t-1 */
double measure_plaquette_sf_weights_improved_boundary_t_minus_1 (int t, double ct) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;

  int x0, x1, x2, x3;

  x0 = t-1;
  
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	
	for (mu1=0;mu1<3;mu1++){
	  
	  ix1=g_iup[ix][mu1];
	  
	  for (mu2=mu1+1;mu2<4;mu2++){ 
	    
	    ix2=g_iup[ix][mu2];
	    
	    v=&g_gauge_field[ix][mu1];
	    w=&g_gauge_field[ix1][mu2];
	    
	    _su3_times_su3(pr1,*v,*w);
	    
	    
	    v=&g_gauge_field[ix][mu2];
	    w=&g_gauge_field[ix2][mu1];
	    
	    _su3_times_su3(pr2,*v,*w);
	    
	    
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    
	    if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	      
	      ac *= ct;
	      
	    }
	    
	    sum += ac;
	    
	  }
	}
	
      }
    }
  }
  
  ga = sum*2.0;

  return ga;
  
}

/*---------------------------------------------------------------------------------------------------------------*/


/* the next functions are used to calculate the Iwasaki action
   it is done "hard-coded"
   but all the kinds of plaquettes and rectangles have a corresponding coefficient
   which is passed as a parameter of the functions so we can change them automatically */

/* plaquette for Iwasaki */
double measure_plaquette_sf_iwasaki(int t, double cs, double ct, double c0) {

  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac; 
  double sum = 0;


  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<3;mu1++){
      
      ix1=g_iup[ix][mu1];
      
      for (mu2=mu1+1;mu2<4;mu2++){ 
	
	ix2=g_iup[ix][mu2];
	
	v=&g_gauge_field[ix][mu1];
	w=&g_gauge_field[ix1][mu2];
	
	_su3_times_su3(pr1,*v,*w);
	
	
	v=&g_gauge_field[ix][mu2];
	w=&g_gauge_field[ix2][mu1];
	
	_su3_times_su3(pr2,*v,*w);
	
	
	_trace_su3_times_su3d(ac,pr1,pr2);

	if (g_t[ix] == 0) {
	  
	  if (mu1 != 0 && mu2 != 0) {

	    ac *= cs;

	  }

	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {

	    ac *= ct;

	  }
	
	} 
	
	
	else if ((g_t[ix] == t) && (mu1 != 0 && mu2 != 0)) {

	    ac *= cs;
	    
	}

	else if ((g_t[ix] == (t-1)) && ((mu1 == 0 || mu2 == 0) && mu1 != mu2)) {

	    ac *= ct;
	    
	}

	else {

	  ac *= c0;

	}
	
	sum += ac;
	
      }
    }
    
  }
  
  ga = sum*2.0;

  return ga;
  
}

/* rectangle for Iwasaki */
double measure_rectangle_sf_iwasaki(int t, double c1, double c1_ss, double c1_tss, double c1_tts) {

  int ix,ix1,ix2,mu1,mu2;
  int ix12,ix22;                           /* My changes */
  static su3 pr_r1,pr_r2,pr_r11,pr_r22;    /* My changes */
  static double ga,ac; 
  su3 *v_r1,*w_r1,*z_r1;                   /* My changes */
  su3 *v_r2,*w_r2,*z_r2;                   /* My changes */
  double sum = 0;                          /* My changes */


  for (ix=0;ix<VOLUME;ix++){
    
    for (mu1=0;mu1<4;mu1++){
      
      ix1 = g_iup[ix][mu1];
      
      for (mu2=0;mu2<4;mu2++){
	if( mu1 != mu2 ){
	  
	  ix2  = g_iup[ix][mu2];
	  ix12 = g_iup[ix1][mu2];
	  ix22 = g_iup[ix2][mu2];
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu2];
	  z_r1 = &g_gauge_field[ix12][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu2];
	  z_r2 = &g_gauge_field[ix22][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr_r22,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr_r11,pr_r22);


	  if ((g_t[ix] == 0 || g_t[ix] == t) && (mu1 != 0 && mu2 != 0)) {

	    ac *= c1_ss;
	  
	  }

	  else if (g_t[ix] == 0 && mu1 == 0) {/* 1 movement in t <=> 2 links on a (time) boundary */

	    ac *= c1_tss;
	  
	  }

	  else if (g_t[ix] == 0 && mu2 == 0) {/* 2 movement in t <=> 1 links on a (time) boundary */

	    ac *= c1_tts;
	  
	  }

	  else if (g_t[ix] == (t-1) && mu1 == 0) {/* 1 movement in t <=> 2 links on a (time) boundary */

	    ac *= c1_tss;
	  
	  }

	  else if (g_t[ix] == (t-2) && mu2 == 0) {/* 2 movement in t <=> 1 links on a (time) boundary */

	    ac *= c1_tts;
	  
	  }

	  else {

	    ac *= c1;

	  }
	  
	  sum += ac;
	}
      }
    }
  }
  
  ga = sum*2.0;

  return ga;

}



/***** ACTIONS *****/

/*----------------------------------------------------------*/

/*** for PBC ***/

/* standard Wilson */
double measure_wilson_action(double beta) {

  double plaquette;
  double wilson;

  plaquette = measure_plaquette();

  wilson = - (beta/(2.*3.)) * plaquette;
  //wilson = beta * (6.*VOLUME*g_nproc - plaquette);

  return wilson;
}

/* Iwasaki */
double measure_iwasaki_action(double beta, double c0, double c1) {

  double plaquette;
  double rectangle;
  double iwasaki;

  plaquette = measure_plaquette();

  rectangle = measure_rectangle();

  iwasaki = - (beta/(2.*3.)) * ( c0*plaquette + c1*rectangle );

  return iwasaki;
}

/*----------------------------------------------------------*/


/*** SF boundary conditions ***/

/* standard Wilson action for SF b.c. without improvement coefficients "hard-coded" */
double measure_wilson_action_sf(int t, double beta) {

  double plaquette;
  double wilson;

  plaquette = measure_plaquette_sf_weights(t);

  wilson = - (beta/(2.*3.)) * plaquette;

  return wilson;
}

/* standard Wilson action for SF b.c. WITH improvement coefficients "hard-coded" */
double measure_wilson_action_sf_weights_improvement(int t, double beta, double cs, double ct) {

  double plaquette;
  double wilson;

  plaquette = measure_plaquette_sf_weights_improvement(t, cs, ct);

  wilson = - (beta/(2.*3.)) * plaquette;

  return wilson;
}


/* it is the same as "measure_wilson_action_sf":
standard Wilson action for SF b.c. without improvement coefficients
but "not hard-coded" */
double measure_wilson_action_sf_separate_boundary(int t, double beta) {

  double plaquette;
  double wilson;

  plaquette = measure_plaquette_sf_weights_bulk(t)
    + measure_plaquette_sf_weights_boundary_0() + measure_plaquette_sf_weights_boundary_t(t);

  wilson = - (beta/(2.*3.)) * plaquette;

  return wilson;
}

/* it is the same as "measure_wilson_action_sf_weights_improvement":
standard Wilson action for SF b.c. WITH (alpha collab.) improvement coefficients
but "not hard-coded" */
double measure_wilson_action_sf_weights_improvement_separate_boundary(int t, double beta, double cs, double ct) {

  double plaquette;
  double wilson;

  plaquette = measure_plaquette_sf_weights_improved_bulk(t) +
    measure_plaquette_sf_weights_improved_boundary_0(cs, ct) +
    measure_plaquette_sf_weights_improved_boundary_t(t, cs) +
    measure_plaquette_sf_weights_improved_boundary_t_minus_1(t, ct);

  wilson = - (beta/(2.*3.)) * plaquette;

  return wilson;
}

/* Iwasaki action with SF b.c. "hard-coded" */
double measure_iwasaki_action_sf(int t, double beta, double cs, double ct, double c0,
				   double c1, double c1_ss, double c1_tss, double c1_tts) {

  double plaquette;
  double rectangle;
  double iwasaki;

  plaquette = measure_plaquette_sf_iwasaki(t, cs, ct, c0);

  rectangle = measure_rectangle_sf_iwasaki(t, c1, c1_ss, c1_tss, c1_tts);

  if (c0 == 1.) {
    iwasaki = - (beta/(2.*3.)) * plaquette;
  }
  else if (c0 == 0.) {
    iwasaki = - (beta/(2.*3.)) * rectangle;
  }
  else {
    iwasaki = - (beta/(2.*3.)) * ( plaquette + rectangle );
  }
  
  return iwasaki;
}




/*** FUNCTIONS NEEDED FOR THE BACKGROUND FIELD ACTION and
     BACKGROUND FIELD ACTION AND DERIVATIVE WITH RESPECT TO ETA ***/


/* it calculates an su3 matrix "u" which is gonna be the partial with respect to eta of the
   (lattice) spatially constant abelian field "C_k"*/
#define _su3_partial_eta_spatially_constant_abelian_field_phi(u)	\
  (u).c00.re = cos(1./LX);						\
  (u).c00.im = sin(1./LX);						\
  (u).c01.re = 0.0;							\
  (u).c01.im = 0.0;							\
  (u).c02.re = 0.0;							\
  (u).c02.im = 0.0;							\
  (u).c10.re = 0.0;							\
  (u).c10.im = 0.0;							\
  (u).c11.re = cos((-1./2.)/LX);					\
  (u).c11.im = sin((-1./2.)/LX);					\
  (u).c12.re = 0.0;							\
  (u).c12.im = 0.0;							\
  (u).c20.re = 0.0;							\
  (u).c20.im = 0.0;							\
  (u).c21.re = 0.0;							\
  (u).c21.im = 0.0;							\
  (u).c22.re = cos((-1./2.)/LX);					\
  (u).c22.im = sin((-1./2.)/LX);


/* it calculates an su3 matrix "u" which is gonna be the partial with respect to eta of the
   (lattice) spatially constant abelian field "(C_k)^prime"*/
#define _su3_partial_eta_spatially_constant_abelian_field_phi_prime(u)	\
  (u).c00.re = cos(-1./LX);						\
  (u).c00.im = sin(-1./LX);						\
  (u).c01.re = 0.0;							\
  (u).c01.im = 0.0;							\
  (u).c02.re = 0.0;							\
  (u).c02.im = 0.0;							\
  (u).c10.re = 0.0;							\
  (u).c10.im = 0.0;							\
  (u).c11.re = cos((1./2.)/LX);						\
  (u).c11.im = sin((1./2.)/LX);						\
  (u).c12.re = 0.0;							\
  (u).c12.im = 0.0;							\
  (u).c20.re = 0.0;							\
  (u).c20.im = 0.0;							\
  (u).c21.re = 0.0;							\
  (u).c21.im = 0.0;							\
  (u).c22.re = cos((1./2.)/LX);						\
  (u).c22.im = sin((1./2.)/LX);


/*--------------------------------------------------------------------------------------------------*/

/** PLAQUETTE (only) **/

/* this function defines the (continuum) constant abelian induced background field "B_{mu}(x)" */
/* that is, the minimal action configuration when SF b.c. are considered
 (note that in an infinite extent, p.b.c., the minimal action configuration is A_{mu}(x)=0) */
void induced_continuum_background(su3 **b, int t, double eta) {
  
  int ix;
  double pi;
  double phi1_0, phi2_0, phi3_0;
  double phi1_T, phi2_T, phi3_T;
  double p1, p2, p3;

  pi = acos(-1.);

  phi1_0 = eta - pi/3.0;
  phi2_0 = - 0.5 * eta;
  phi3_0 = - 0.5 * eta + pi/3.0;
  
  phi1_T = - phi1_0 - (4.0*pi)/3.0;
  phi2_T = - phi3_0 + (2.0*pi)/3.0;
  phi3_T = - phi2_0 + (2.0*pi)/3.0;

  phi1_0 /= (double)LX; 
  phi2_0 /= (double)LX; 
  phi3_0 /= (double)LX; 
  
  phi1_T /= (double)LX; 
  phi2_T /= (double)LX; 
  phi3_T /= (double)LX;

  phi1_0 /= (double)t; 
  phi2_0 /= (double)t; 
  phi3_0 /= (double)t; 
  
  phi1_T /= (double)t; 
  phi2_T /= (double)t; 
  phi3_T /= (double)t; 

  
  for (ix=0;ix<VOLUME;ix++){

    p1 = g_t[ix]*phi1_T + ((double)t - g_t[ix])*phi1_0;
    p2 = g_t[ix]*phi2_T + ((double)t - g_t[ix])*phi2_0;
    p3 = g_t[ix]*phi3_T + ((double)t - g_t[ix])*phi3_0;

    _su3_zero(b[ix][0]);
    _su3_spatially_constant_abelian_field_continuum(b[ix][1], p1, p2, p3);
    _su3_spatially_constant_abelian_field_continuum(b[ix][2], p1, p2, p3);
    _su3_spatially_constant_abelian_field_continuum(b[ix][3], p1, p2, p3);


  }
}


/* this function defines the lattice constant abelian induced background field:
   "V(x,mu)=exp{aB_{mu}(x)}" of the plaquette action
   (note that with pbc this would be U(x,mu)=Identity) */

/* Note that this configuration will be also the minimal one in the case of the Iwasaki action
   when option B of the paper hep-lat/9808007 is chosen */

/* For other choice of the weight factors (not option B), this will NOT be the minimal action configuration.
   Actually, besides in this case B, the analytical expression will not be know so,
   only a numerical determination of the minimal configuration will be possible */
void induced_lattice_background(su3 **v, int t, double eta) { 

  int ix;
  double pi;
  double phi1_0, phi2_0, phi3_0;
  double phi1_T, phi2_T, phi3_T;
  double p1, p2, p3;

  pi = acos(-1.);

  phi1_0 = eta - pi/3.0;
  phi2_0 = - 0.5 * eta;
  phi3_0 = - 0.5 * eta + pi/3.0;
  
  phi1_T = - phi1_0 - (4.0*pi)/3.0;
  phi2_T = - phi3_0 + (2.0*pi)/3.0;
  phi3_T = - phi2_0 + (2.0*pi)/3.0;

  phi1_0 /= (double)LX; 
  phi2_0 /= (double)LX; 
  phi3_0 /= (double)LX; 
  
  phi1_T /= (double)LX; 
  phi2_T /= (double)LX; 
  phi3_T /= (double)LX;

  phi1_0 /= (double)t; 
  phi2_0 /= (double)t; 
  phi3_0 /= (double)t; 
  
  phi1_T /= (double)t; 
  phi2_T /= (double)t; 
  phi3_T /= (double)t;



  for (ix=0;ix<VOLUME;ix++){

    p1 = g_t[ix]*phi1_T + ((double)t - g_t[ix])*phi1_0;
    p2 = g_t[ix]*phi2_T + ((double)t - g_t[ix])*phi2_0;
    p3 = g_t[ix]*phi3_T + ((double)t - g_t[ix])*phi3_0;


    if(g_t[ix] == t) {
      
      _su3_zero(v[ix][0]);
      
    }
    
    else {
      
      _su3_one(v[ix][0]);
      
    }
    
    _su3_spatially_constant_abelian_field(v[ix][1], p1, p2, p3);
    _su3_spatially_constant_abelian_field(v[ix][2], p1, p2, p3);
    _su3_spatially_constant_abelian_field(v[ix][3], p1, p2, p3);
   
  }
}




/* the next function gives us the classical lattice background field (V) "plaquette action S[V]" with SF b.c. */
/* implementation of the analytical expression taken from Rainer's notes in Schladming (eq.71)*/
double lattice_background_plaquette_action_sf(int t, double beta, double ct, double eta) {

  double pi;
  double phi1_0, phi2_0, phi3_0;
  double phi1_T, phi2_T, phi3_T;
  double plaquette_back;
  double factor1, factor2, factor3, factor;
  double factor_a, factor_b, factor_c;

  pi = acos(-1.);

  phi1_0 = eta - pi/3.;
  phi2_0 = - 0.5 * eta;
  phi3_0 = - 0.5 * eta + pi/3.;
  
  phi1_T = - phi1_0 - (4.*pi)/3.;
  phi2_T = - phi3_0 + (2.*pi)/3.;
  phi3_T = - phi2_0 + (2.*pi)/3.;
  

  factor1 = (1. - (1. - ct)*(2./(double)t));
  factor2 = (beta*(double)LX*(double)LX*(double)LX*(double)t)/2.;

  factor3 = 2.*2.*factor1*factor2;

  factor = 1./(2.*(double)LX*(double)t);

  factor_a = factor*(phi1_T - phi1_0);
  factor_b = factor*(phi2_T - phi2_0);
  factor_c = factor*(phi3_T - phi3_0);

  plaquette_back = factor3*(sin(factor_a)*sin(factor_a) + sin(factor_b)*sin(factor_b) + sin(factor_c)*sin(factor_c));
  
  
  return plaquette_back;
}

/* the next function gives us the "leading order" of the "effective action" with SF b.c. */
double lattice_lo_effective_plaquette_action_sf(int t, double beta, double ct, double eta) {

  double eff_lo;
  double plaq_back;

  plaq_back = lattice_background_plaquette_action_sf(t, beta, ct, eta);

  eff_lo = (6./beta)*plaq_back;

  return eff_lo;
}

/* the next function gives us the "partial derivative" with respect to "eta"
   of the classical lattice background field (V) "plaquette action S[V]" with SF b.c. */
/* implementation of the derivative of the analytical expression taken from Rainer's notes in Schladming (eq.71)*/
double partial_lattice_background_plaquette_action_sf(int t, double beta, double ct, double eta) {

  double pi;
  double phi1_0, phi2_0, phi3_0;
  double phi1_T, phi2_T, phi3_T;
  double phi1_0p, phi2_0p, phi3_0p;
  double phi1_Tp, phi2_Tp, phi3_Tp;
  double partial_plaquette_back;
  double factor1, factor2, factor3, factor;
  double factor_a, factor_b, factor_c;
  double diff_a, diff_b, diff_c;
  double a,b,c;

  pi = acos(-1.);


  phi1_0 = eta - pi/3.0;
  phi2_0 = - 0.5 * eta;
  phi3_0 = - 0.5 * eta + pi/3.0;
  
  phi1_T = - phi1_0 - (4.0*pi)/3.0;
  phi2_T = - phi3_0 + (2.0*pi)/3.0;
  phi3_T = - phi2_0 + (2.0*pi)/3.0;


  phi1_0p = 1.0;
  phi2_0p = - 0.5;
  phi3_0p = - 0.5;
  
  phi1_Tp = - phi1_0p;
  phi2_Tp = - phi3_0p;
  phi3_Tp = - phi2_0p;
  


  factor1 = (1. - (1. - ct)*(2./(double)t));

  printf("factor1 = %e \n", factor1);

  factor2 = 2.*beta*(double)LX*(double)LX;

  factor3 = factor1*factor2;


  factor = 1./(2.*(double)LX*(double)t);

  factor_a = factor*(phi1_T - phi1_0);
  factor_b = factor*(phi2_T - phi2_0);
  factor_c = factor*(phi3_T - phi3_0);

  diff_a = phi1_Tp - phi1_0p;
  diff_b = phi2_Tp - phi2_0p;
  diff_c = phi3_Tp - phi3_0p;

  a = sin(factor_a)*cos(factor_a)*diff_a;
  b = sin(factor_b)*cos(factor_b)*diff_b;
  c = sin(factor_c)*cos(factor_c)*diff_c;

  partial_plaquette_back = factor3*(a + b + c);
  
  
  return partial_plaquette_back;
}


/* the next function gives us the "derivative" with respect to "eta"
   of the "leading order" of the "effective action" with SF b.c. */
double partial_lattice_lo_effective_plaquette_action_sf(int t, double beta, double ct, double eta) {

  double partial_eff_lo;
  double partial_plaq_back;

  partial_plaq_back = partial_lattice_background_plaquette_action_sf(t, beta, ct, eta);

  partial_eff_lo = (6./beta)*partial_plaq_back;

  return partial_eff_lo;
}


/*--------------------------------------------------------------------------------------------------*/


#if 0
/** IWASAKI **/
/* STILL TO DO!!!!!!!!!! */
/* these are all implementations of analytical expressions */
/* note that we ONLY know analytically a solution for the EOM for the lattice background field
   of the IWASAKI acion in the case we choose Option_B of the paper: hep-lat/9808007, 
   which in this case, actually coincides with the expression obtained by the alpha_collab. for 
   the plaquette Wilson action with sf */
/* Therefore, whenever we are not considering this case B,
   the minimal lattice action (Iwasaki) configuration V must be found only numerically */

double lattice_background_******_action_sf(int t, double beta, double ct, double eta) {

}


double lattice_lo_effective_*****_action_sf(int t, double beta, double ct, double eta) {

}

double partial_lattice_background_****_action_sf(int t, double beta, double ct, double eta) {

}

double partial_lattice_lo_effective_****_action_sf(int t, double beta, double ct, double eta) {

}
#endif

/*-------------------------------------------------------------------------------------------------*/




/*** DEFINITION OF THE RUNNING COUPLING ***/


/* it calculates an su3 matrix "u" which is gonna be the eighth-generator of SU(3) "lambda_8" */
#define _su3_lambda_8(u)					\
  (u).c00.re = 1.0;						\
  (u).c00.im = 0.0;						\
  (u).c01.re = 0.0;						\
  (u).c01.im = 0.0;						\
  (u).c02.re = 0.0;						\
  (u).c02.im = 0.0;						\
  (u).c10.re = 0.0;						\
  (u).c10.im = 0.0;						\
  (u).c11.re =-0.5;						\
  (u).c11.im = 0.0;						\
  (u).c12.re = 0.0;						\
  (u).c12.im = 0.0;						\
  (u).c20.re = 0.0;						\
  (u).c20.im = 0.0;						\
  (u).c21.re = 0.0;						\
  (u).c21.im = 0.0;						\
  (u).c22.re =-0.5;						\
  (u).c22.im = 0.0;


/* it calculates an su3 matrix "u" which is gonna be the eighth-generator of SU(3) "lambda_8" times "i" */
#define _su3_i_times_lambda_8(u)				\
  (u).c00.re = 0.0;						\
  (u).c00.im = 1.0;						\
  (u).c01.re = 0.0;						\
  (u).c01.im = 0.0;						\
  (u).c02.re = 0.0;						\
  (u).c02.im = 0.0;						\
  (u).c10.re = 0.0;						\
  (u).c10.im = 0.0;						\
  (u).c11.re = 0.0;						\
  (u).c11.im =-0.5;						\
  (u).c12.re = 0.0;						\
  (u).c12.im = 0.0;						\
  (u).c20.re = 0.0;						\
  (u).c20.im = 0.0;						\
  (u).c21.re = 0.0;						\
  (u).c21.im = 0.0;						\
  (u).c22.re = 0.0;						\
  (u).c22.im =-0.5;



/*------------------------------------------------------------------------------------------------------------*/
 
/*
void testfunc() {
  static  su3 a,b,c;
  _su3_times_su3(a,b,c);
}
*/



/* it gives the expression of the "derivative" of the "PLAQUETTE" with SF b.c.
   with respect to the background field parameter "eta" */
/* it has been taken from Rainer's notes in Schladming (eq.73) (we've checked and gotten the same formula) */
/* WARNING: this function is only valid if we are considering U!=V */
double partial_plaquette_sf_respect_to_eta(int t, double ct) {
  
  
  int ix,ix1,ix2,mu1,mu2;  
  static su3 pr1,pr2; 
  su3 *v,*w;
  static double ga,ac;

  int x0, x1, x2, x3; 
  double sum_0 = 0, sum_tminus1 = 0;  
  double ga_0, ga_tminus1, ga_int;
  static su3 pr11;  
  static su3 i_lambda8;
  
  /*printf("%ld\n", (((unsigned long int)(&pr11))%16) );
    printf("%ld\n", (((unsigned long int)(&i_lambda8))%16) );*/
  
  _su3_i_times_lambda_8(i_lambda8);
  
  mu2 = 0;
  
  /* loop to calculate E_{k}^{8}(\vec{x}) */
  x0 = 0;
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=1;mu1<4;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  
	  v=&g_gauge_field[ix][mu1];
	  w=&g_gauge_field[ix1][mu2];
	  
	  _su3_times_su3(pr11,*v,*w);
	  
	  _su3_times_su3(pr1,i_lambda8, pr11);
	  
	  v=&g_gauge_field[ix][mu2];
	  w=&g_gauge_field[ix2][mu1];
	  
	  _su3_times_su3(pr2,*v,*w);
	  
	    _trace_su3_times_su3d(ac,pr1,pr2);
	    
	    sum_0 += ac;
	    
	    
	}
	
      }
    }
  }
  
  ga_0 = sum_0;

  
  /* loop to calculate (E_{k}^{8})^{prime}(\vec{x}) */
  x0 = t-1;
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=1;mu1<4;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  
	  v=&g_gauge_field[ix][mu1];
	  w=&g_gauge_field[ix1][mu2];
	  
	  _su3_times_su3(pr11,*v,*w);
	  _su3_times_su3(pr1,i_lambda8,pr11);
	  
	  v=&g_gauge_field[ix][mu2];
	  w=&g_gauge_field[ix2][mu1];
	  
	  _su3_times_su3(pr2,*v,*w);
	  
	  _trace_su3_times_su3d(ac,pr1,pr2);
	  
	  sum_tminus1 += ac;
	  
	}
	
      }
    }
  }
  
  ga_tminus1 = sum_tminus1;
  
  ga_int = ct*(ga_0 + ga_tminus1);
  
  /*ga_int = ct*ga_0 - ct*ga_tminus1; */ /* this was WRONG!!! */
  ga = ga_int;
  
#if 0
  
    printf("\n"); fflush(stdout);
    printf("sum_0 = %e \n",sum_0); fflush(stdout);
    printf("ga_0 = %e \n",ga_0); fflush(stdout);
    
    printf("\n"); fflush(stdout);
    printf("sum_tminus1 = %e \n",sum_tminus1); fflush(stdout);
    printf("ga_tminus1 = %e \n",ga_tminus1); fflush(stdout);
    
    printf("\n"); fflush(stdout);
    printf("ct = %e \n", ct); fflush(stdout);
    
    printf("\n"); fflush(stdout);  /* if I eliminate this line ===> segmentation fault */
    printf("ga = %e \n", ga); fflush(stdout);
    printf("\n"); fflush(stdout);
#endif

  return ga;
  
}

/*------------------------------------------------------------------------------------------------------------*/


/* it gives the expression of the "derivative" of the "RECTANGLE" with SF b.c.
   with respect to the background field parameter "eta" */
/* it has been derived similarly as has been done for the plquette in the previous case */
double partial_rectangle_sf_respect_to_eta(int t, double c1_tss, double c1_tts) {
  
  int ix,ix1,ix2,mu1,mu2;
  static su3 pr1,pr2; 
  static double ga,ac;
  
  int x0, x1, x2, x3;
  int ix12,ix11,ix22; 
  static su3 i_lambda8;
  su3 *v_r1,*w_r1,*z_r1;
  su3 *v_r2,*w_r2,*z_r2;
  static su3 pr_r1,pr_r2,pr_r11;
  double sum_0_tss=0, sum_0_tts=0, sum_tminus1_tss=0, sum_tminus2_tts=0;  
  double ga_0_tss, ga_0_tts,  ga_tminus1_tss, ga_tminus2_tts;
  
  
  _su3_i_times_lambda_8(i_lambda8);
  
  mu2 = 0;

  
 /* loops to calculate E_{k}^{8}(\vec{x}): "R_tss" and "R_tts" at 0 */ 
  x0 = 0;
  
  
  /* first contribution to E_{k}^{8}(\vec{x}): R_tss */
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=1;mu1<4;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  ix12=g_iup[ix1][mu2];
	  ix11=g_iup[ix1][mu1];	
	  ix22=g_iup[ix2][mu2];	
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu1];
	  z_r1 = &g_gauge_field[ix11][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  _su3_times_su3(pr1,i_lambda8,pr_r11);	    
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu1];
	  z_r2 = &g_gauge_field[ix12][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr2,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr1,pr2);
	  
	  
	  sum_0_tss += ac;
	  
	}
	
      }
    }
  }
  
  ga_0_tss = sum_0_tss;
  
  
  /* second contribution to E_{k}^{8}(\vec{x}): R_tts */
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=1;mu1<4;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  ix12=g_iup[ix1][mu2];
	  ix11=g_iup[ix1][mu1];	
	  ix22=g_iup[ix2][mu2];	
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu2];
	  z_r1 = &g_gauge_field[ix12][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  _su3_times_su3(pr1,i_lambda8,pr_r11);	    
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu2];
	  z_r2 = &g_gauge_field[ix22][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr2,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr1,pr2);
	  
	  
	  sum_0_tts += ac;
	  
	}
	
      }
    }
  }
  
  ga_0_tts = sum_0_tts;



    /* loop to calculate (E_{k}^{8})^{prime}(\vec{x}): "R_tss" and "R_tts" at t */

    /* first contribution to (E_{k}^{8})^{prime}(\vec{x}): R_tss */
  x0 = t-1;
  
  for (x1 = 0 ; x1 < LX ; x1++) {
    
    for (x2 = 0 ; x2 < LY ; x2++) {
      
      for (x3 = 0 ; x3 < LZ ; x3++) {
	
	ix = Index(x0, x1, x2, x3);
	
	for (mu1=1;mu1<4;mu1++) {
	  
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  ix12=g_iup[ix1][mu2];
	  ix11=g_iup[ix1][mu1];	
	  ix22=g_iup[ix2][mu2];	
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu1];
	  z_r1 = &g_gauge_field[ix11][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  _su3_times_su3(pr1,i_lambda8,pr_r11);	    
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu1];
	  z_r2 = &g_gauge_field[ix12][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr2,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr1,pr2);
	  
	  
	  sum_tminus1_tss += ac;
	  
	}
	
      }
    }
  }
  
  ga_tminus1_tss = sum_tminus1_tss;
  

    /* second contribution to (E_{k}^{8})^{prime}(\vec{x}): R_tts */
    x0 = t-2;

    for (x1 = 0 ; x1 < LX ; x1++) {
      
      for (x2 = 0 ; x2 < LY ; x2++) {
	
	for (x3 = 0 ; x3 < LZ ; x3++) {
	  
	  ix = Index(x0, x1, x2, x3);
	  
	  for (mu1=1;mu1<4;mu1++) {
	    
	  ix1=g_iup[ix][mu1];
	  ix2=g_iup[ix][mu2];
	  ix12=g_iup[ix1][mu2];
	  ix11=g_iup[ix1][mu1];	
	  ix22=g_iup[ix2][mu2];	
	  
	  
	  v_r1 = &g_gauge_field[ix][mu1];
	  w_r1 = &g_gauge_field[ix1][mu2];
	  z_r1 = &g_gauge_field[ix12][mu2];
	  
	  _su3_times_su3(pr_r1,*v_r1,*w_r1);
	  _su3_times_su3(pr_r11,pr_r1,*z_r1);
	  _su3_times_su3(pr1,i_lambda8,pr_r11);	    
	  
	  
	  v_r2 = &g_gauge_field[ix][mu2];
	  w_r2 = &g_gauge_field[ix2][mu2];
	  z_r2 = &g_gauge_field[ix22][mu1];
	  
	  _su3_times_su3(pr_r2,*v_r2,*w_r2);
	  _su3_times_su3(pr2,pr_r2,*z_r2);
	  
	  
	  _trace_su3_times_su3d(ac,pr1,pr2);
	  
	  
	  sum_tminus2_tts += ac;

	    
	    }
	  
	}
      }
    }

    ga_tminus2_tts = sum_tminus2_tts;



    /* ga = c1_tss*(ga_0_tss - ga_tminus1_tss) + c1_tts*(ga_0_tts - ga_tminus2_tts); */ /* it was also WRONG!!! */
    ga = 2.*c1_tss*(ga_0_tss + ga_tminus1_tss) + c1_tts*(ga_0_tts + ga_tminus2_tts);


  return ga;  
  
  
}

/*------------------------------------------------------------------------------------------------------------*/

/* the next function gives us the "derivative" of the "IWASAKI action" with SF b.c.
   with respect to the background field parameter "eta" */
/* WARNING: this function is only valid if we are considering U!=V */
double partial_iwasaki_action_sf_respect_to_eta(int t, double beta, double cs, double ct, double c0,
						double c1, double c1_ss, double c1_tss, double c1_tts) {

  double partial_plaquette;
  double partial_rectangle;
  double partial_iwasaki;


  partial_plaquette = partial_plaquette_sf_respect_to_eta(t, ct);
  partial_rectangle = partial_rectangle_sf_respect_to_eta(t, c1_tss, c1_tts);
  
  if (c0 == 1.) {
    partial_iwasaki = - (beta/(3.*(double)LX)) * partial_plaquette;
  }
  else if (c0 == 0.) {
    partial_iwasaki = - (beta/(3.*(double)LX)) * partial_rectangle;
  }
  else {
    partial_iwasaki = - (beta/(3.*(double)LX)) * ( partial_plaquette + partial_rectangle );
  }
  
  return partial_iwasaki;
}


/********************************************************************************************************/
/********************************************************************************************************/

#define _su3_nan(u) \
   (u).c00.re=1./0.0; \
   (u).c00.im=1./0.0; \
   (u).c01.re=1./0.0; \
   (u).c01.im=1./0.0; \
   (u).c02.re=1./0.0; \
   (u).c02.im=1./0.0; \
   (u).c10.re=1./0.0; \
   (u).c10.im=1./0.0; \
   (u).c11.re=1./0.0; \
   (u).c11.im=1./0.0; \
   (u).c12.re=1./0.0; \
   (u).c12.im=1./0.0; \
   (u).c20.re=1./0.0; \
   (u).c20.im=1./0.0; \
   (u).c21.re=1./0.0; \
   (u).c21.im=1./0.0; \
   (u).c22.re=1./0.0; \
   (u).c22.im=1./0.0;


/* the next function sets
all the gauge links in the time direction (from t on) to nan.
Note that the rest of the links at the boundaries (spatial links) are not yet touched here */
void nan_dirichlet_boundary_conditions(int t) {
  
  int ix;
  
  for (ix=0;ix<VOLUME;ix++){
 
    if (g_t[ix] == t) {

      _su3_nan(g_gauge_field[ix][0]);

    }   
  }
}









/********************************************************************************************************/
/********************************************************************************************************/











/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/


#if 0
  double plaquette_energy;
  double rectangle_energy;
  double wilson_action;
  double wilson_action_sepbound;
  double iwasaki_action;
  double partial_iwa;
  double partial_iwasaki_action;



  /* some parameters (Jen) */

  t_bsf = T-1; /* it sets at which time slice I want to put the SF b.c. (end point) --- T = lattice time extent set by Carsten */
  printf("\n"); fflush(stdout);
  printf("t_bsf = %i \n", g_Tbsf); fflush(stdout);

  /* Action "A" of "hep-lat/9808007" */
#if 0
  c_1_ss = 0.0;
  c_1_tss = c_1;
  c_1_tts = c_1;
#endif

  /* Action "B" of "hep-lat/9808007" */
#if 0
  c_1_ss = 0.0;
  c_1_tss = (3./2.)*c_1;
  c_1_tts = c_1;
#endif


  /*** PERIODIC b.c. ***/

  /*compute the energy of the gauge field*/
  /*plaquette_energy = measure_gauge_action();*/
  plaquette_energy = measure_plaquette();
  rectangle_energy = measure_rectangle();
  wilson_action = measure_wilson_action(g_beta);
  iwasaki_action = measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("Periodic boundary conditions: \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);

  }

  /*** PREPARING for SF b.c. ***/

  /* fix Dirichlet bc in time direction */
  dirichlet_boundary_conditions(g_Tbsf);

  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette();
  rectangle_energy = measure_rectangle();
  wilson_action = measure_wilson_action(g_beta);
  iwasaki_action = measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("Dirichlet boundary conditions: \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);

  }

#if 0
  /* fix Dirichlet bc in the time direction and with spatial links to one */
  dirichlet_boundary_conditions_spatial_links_to_one(g_Tbsf);

  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette();
  rectangle_energy = measure_rectangle();
  wilson_action = measure_wilson_action(g_beta);
  iwasaki_action = measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("Dirichlet boundary conditions with spatial links to one: \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);

  }
#endif



  /*** SF BOUNDARY CONDITIONS ***/

  /* sf b.c. abelian field */
  sf_boundary_conditions_spatially_constant_abelian_field(g_Tbsf, g_eta);


#if 0
  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette();
  rectangle_energy = measure_rectangle();
  wilson_action = measure_wilson_action(g_beta);
  iwasaki_action = measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian: \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);
 
 }
#endif

#if 0
  /* sf b.c. abelian field and standard sf weight factors included (only plaquette here) */

  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette_sf_weights(g_Tbsf);
  wilson_action = measure_wilson_action_sf(g_Tbsf, g_beta);
  wilson_action_sepbound = measure_wilson_action_sf_separate_boundary(g_Tbsf, g_beta);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and standard sf weight factors included (only plaquette): \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Wilson action value sep bound is %e\n", wilson_action_sepbound); fflush(stdout);
 
 }
#endif

  /* sf b.c. abelian field and weight factors for O(a)-improvement included (only plaquette here) */

  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette_sf_weights_improvement(g_Tbsf, g_Cs, g_Ct) ;
  wilson_action = measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct);
  wilson_action_sepbound = measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and weight factors for O(a)-improvement included (only plaquette): \n"); fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Wilson action value is %e\n", wilson_action); fflush(stdout);
    printf("The Wilson action value sep bound is %e\n", wilson_action_sepbound); fflush(stdout);
 
 }

  /* sf b.c. abelian field and weight factors for O(a)-improvement included (plaquette and rectangle) */

  /*compute the energy of the gauge field*/
  plaquette_energy = measure_plaquette_sf_iwasaki(g_Tbsf, g_Cs, g_Ct, g_rgi_C0) ;
  rectangle_energy = measure_rectangle_sf_iwasaki(g_Tbsf, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);

  if(g_proc_id == 0) {
    printf("\n"); fflush(stdout);
    printf("SF b.c. abelian and weight factors for O(a)-improvement included (Iwasaki = plaquette and rectangle): \n");
    fflush(stdout);
    printf("The plaquette value is %e\n", plaquette_energy/(3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The rectangle value is %e\n", rectangle_energy/(2.*3.*6.*VOLUME*g_nproc)); fflush(stdout);
    printf("The Iwasaki action value is %e\n", iwasaki_action); fflush(stdout);

 }



  /*** CHECKS: here we calculate "S[V], Gamma[V], S'[V] and Gamma'[V]" in two ways and both should agree ***/
  printf("\n"); fflush(stdout);
  printf("CHECKS: \n"); fflush(stdout);

#if 1
  /* (0): here we have not yet assigned: U = V forall x_0 ==> it should not agree with the (1) and (2) */

  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action = partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_eta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf("\n"); fflush(stdout);
  printf(" Before assigning U=V but SF b.c. \n"); fflush(stdout);
  printf("S[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout); 
  printf("\n"); fflush(stdout);
#endif
  
  
  /* (1): identifying the gauge fields "g_gauge_fields = V" and then calculating the plaquette as usually */
  
  induced_lattice_background(g_gauge_field, g_Tbsf, g_eta);
  
  wilson_action = measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct);
  wilson_action_sepbound = measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct);
  iwasaki_action = measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action = partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf(" Assigning U=V with the functions defined for that and then calculating S[V] from the same functions to calculate the actions as in previous cases \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", wilson_action_sepbound); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", wilson_action ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("G[V] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("S'[V] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("G'[V] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);


  /* obtaine normalization factor by calculation Wilson action for U=1 in all the lattice
   and substract it to the previous result for the action.
  Therefore, it should agree with the result obtained from the analytical expression implemented below */
  set_all_links_to_one_with_dirichlet(g_Tbsf);

  iwasaki_action -= measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  partial_iwasaki_action -= partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts);
  
  printf("\n"); fflush(stdout);
  printf(" Previous case but substracting the normalization factor to the action: \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("Norm - S_sf_iwasaki_notsepb[U,W',W] = %e \n", iwasaki_action); fflush(stdout);
  printf("Norm - G[V] = %e \n", (6./g_beta)*iwasaki_action); fflush(stdout);
  printf("Norm' - S'[V] = %e \n", partial_iwasaki_action); fflush(stdout);
  printf("Norm' - G'[V] = %e \n", (6./g_beta)*partial_iwasaki_action);fflush(stdout);
  printf("\n"); fflush(stdout);


  /* (2): directly from the analytical expression which has been implemente in: */
  printf("\n"); fflush(stdout);
  printf(" Assigning U=V: but directly using the analytical expression of the action S[V] \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  printf("S[V]_analy = %e \n", lattice_background_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("G[V]_analy = %e \n", lattice_lo_effective_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("S'[V]_analy = %e \n", partial_lattice_background_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("G'[V]_analy = %e \n", partial_lattice_lo_effective_plaquette_action_sf(g_Tbsf, g_beta, g_Ct, g_eta)); fflush(stdout);
  printf("\n"); fflush(stdout);


#if 1
  /* obtaine normalization factor by calculation Wilson action for U=1 in all the lattice */
  set_all_links_to_one_with_dirichlet(g_Tbsf);

  printf("\n"); fflush(stdout);
  printf(" Setting U=Id and Dirichlet at x0= 0, t \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  /* The next three prints give me the same result, from 3 different functions.
   The first two functions were cross-checked with Dru ==> they should be right.
   Hoever, the result here obtained still differs to what we obtain by doing the
   differenct between our result (for U=V) and the analytical expression */
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct)); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct) ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts));fflush(stdout); 
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);
#endif



#if 0
  /* obtaine normalization factor by calculation Wilson action for U=1 in all the lattice */
  set_all_links_to_one();

  printf("\n"); fflush(stdout);
  printf(" Setting U=Id \n"); fflush(stdout);
  printf("\n"); fflush(stdout);
  /* For the first case below, pbc, I've gotten the number I expected: "(Nc*12*L^4)/g02".
   Thus, since the function "measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1))" was crosschecked bf with Dru it should be right.
  It somehow tells me that also the function which assigns the gauge fields to one "set_all_links_to_one()" should be right.*/
  printf("S_pbc[U,W',W] = %e \n", measure_iwasaki_action(g_beta, g_rgi_C0, g_rgi_C1)); fflush(stdout);
  /* The next three prints give me the same result, from 3 different functions.
   The first two functions were cross-checked with Dru ==> they should be right.
   Hoever, the result here obtained still differs to what we obtain by doing the
   differenct between our result (for U=V) and the analytical expression */
  printf("S_sf_wilson_sepbound[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement_separate_boundary(g_Tbsf, g_beta, g_Cs, g_Ct)); fflush(stdout);
  printf("S_sf_wilson_notsepbd[U,W',W] = %e \n", measure_wilson_action_sf_weights_improvement(g_Tbsf, g_beta, g_Cs, g_Ct) ); fflush(stdout);
  printf("S_sf_iwasaki_notsepb[U,W',W] = %e \n", measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G[U,W',W] = %e \n", (6./g_beta)*measure_iwasaki_action_sf(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("S'[U,W',W] = %e \n", partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts)); fflush(stdout);
  printf("G'[U,W',W] = %e \n", (6./g_beta)*partial_iwasaki_action_sf_respect_to_eta(g_Tbsf, g_beta, g_Cs, g_Ct, g_rgi_C0, g_rgi_C1, g_C1ss, g_C1tss, g_C1tts));fflush(stdout); 
  printf("\n"); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_bulk = %e \n", measure_plaquette_sf_weights_improved_bulk(g_Tbsf)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_0(cs,ct) = %e \n", measure_plaquette_sf_weights_improved_boundary_0(g_Cs, g_Ct)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t(cs) = %e \n", measure_plaquette_sf_weights_improved_boundary_t(g_Tbsf, g_Cs)); fflush(stdout);
  printf("measure_plaquette_sf_weights_improved_boundary_t_minus_1(ct) = %e \n",  measure_plaquette_sf_weights_improved_boundary_t_minus_1(g_Tbsf, g_Ct)); fflush(stdout);
  printf("\n"); fflush(stdout);
#endif


#endif
