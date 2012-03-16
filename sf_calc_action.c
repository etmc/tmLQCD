/*******************************************
*
* FILE: sf_calc_action.c
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

int Index(const int x0, const int x1, const int x2, const int x3);

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
  (u).c00 = cexp(p1);						\
  (u).c01 = 0.0;						\
  (u).c02 = 0.0;						\
  (u).c10 = 0.0;						\
  (u).c11 = cexp(p2);						\
  (u).c12 = 0.0;						\
  (u).c20 = 0.0;						\
  (u).c21 = 0.0;						\
  (u).c22 = cexp(p3);


/* it calculates an su3 matrix "u" which is gonna be the (continuum) spatially constant abelian field */
#define _su3_spatially_constant_abelian_field_continuum(u, p1, p2, p3)	\
  (u).c00 = p1 * I;							\
  (u).c01 = 0.0;							\
  (u).c02 = 0.0;							\
  (u).c10 = 0.0;							\
  (u).c11 = p2 * I;							\
  (u).c12 = 0.0;							\
  (u).c20 = 0.0;							\
  (u).c21 = 0.0;							\
  (u).c22 = p3 * I;


/* it just prints on the screen the su3 matrix "u" */
void print_su3_matrix (su3 u) {
  
  printf(" %e i %e  %e i %e  %e i %e \n",
	 creal(u.c00),  cimag(u.c00),  creal(u.c01),  cimag(u.c01),  creal(u.c02),  cimag(u.c02));
  printf(" %e i %e  %e i %e  %e i %e \n",
	 creal(u.c10),  cimag(u.c10),  creal(u.c11),  cimag(u.c11),  creal(u.c12),  cimag(u.c12));
  printf(" %e i %e  %e i %e  %e i %e \n",
	 creal(u.c20),  cimag(u.c20),  creal(u.c21),  cimag(u.c21),  creal(u.c22),  cimag(u.c22));
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
	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	    
	    ac = 0.;
	    
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
	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	    
	  ac = 0.;
	    
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
	    if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	      
	      ac = 0.;
	      
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
	    if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	      
	      ac = 0.;
	      
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
	else if (g_t[ix] == t) {
	  
	  if (mu1 != 0 && mu2 != 0) {
	    ac *= cs;
	  }
	  if ((mu1 == 0 || mu2 == 0) && mu1 != mu2) {
	    ac = 0.;
	  }
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
	  else if (g_t[ix] == (t-1) && mu2 == 0) {/* 2 movement in t <=> 1 links on a (time) boundary */

	    ac = 0.;
	  
	  }

	  else if (g_t[ix] == (t-2) && mu2 == 0) {/* 2 movement in t <=> 1 links on a (time) boundary */

	    ac *= c1_tts;
	  
	  }
	  else if (g_t[ix] == t && (mu1 == 0 || mu2 == 0)) {/* out of the lattice on the right side */

	    ac = 0.;
	  
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
  /* wilson = beta * (6.*VOLUME*g_nproc - plaquette); */

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


  iwasaki = - (beta/(2.*3.)) * ( plaquette + rectangle );  
  
  return iwasaki;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/


/*** FUNCTIONS NEEDED FOR THE BACKGROUND FIELD ACTION and DERIVATIVE WITH RESPECT TO ETA ***/


/* it calculates an su3 matrix "u" which is gonna be the partial with respect to eta of the
   (lattice) spatially constant abelian field "C_k"*/
#define _su3_partial_eta_spatially_constant_abelian_field_phi(u)	\
  (u).c00 = cexp(1./LX);						\
  (u).c01 = 0.0;							\
  (u).c02 = 0.0;							\
  (u).c10 = 0.0;							\
  (u).c11 = cexp((-1./2.)/LX);						\
  (u).c12 = 0.0;							\
  (u).c20 = 0.0;							\
  (u).c21 = 0.0;							\
  (u).c22 = cexp((-1./2.)/LX);						\


/* it calculates an su3 matrix "u" which is gonna be the partial with respect to eta of the
   (lattice) spatially constant abelian field "(C_k)^prime"*/
#define _su3_partial_eta_spatially_constant_abelian_field_phi_prime(u)	\
  (u).c00 = cexp(-1./LX);						\
  (u).c01 = 0.0;							\
  (u).c02 = 0.0;							\
  (u).c10 = 0.0;							\
  (u).c11 = cexp((1./2.)/LX);						\
  (u).c12 = 0.0;							\
  (u).c20 = 0.0;							\
  (u).c21 = 0.0;							\
  (u).c22 = cexp((1./2.)/LX);						\


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
/* NOTE: there was a sign mistake in Rainer's notes!!! */
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
/* NOTE: there was a sign mistake in Rainer's notes!!! */
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

/* the next function gives the normalisation factor "K" for the definition of the coupling constant
   when considering only the plaquette, as taken from the ALPHA SU(3) paper */
double partial_lattice_lo_effective_plaquette_action_sf_k(int t, double beta, double ct, double eta) {

  double pi;
  double factor1, factor;
  double partial_eff_lo;

  pi = acos(-1.);
  
  factor1 = 12.*(double)LX*(double)LX;
  
  factor = (1./((double)LX*(double)t))*(eta + pi/3.0);

  partial_eff_lo = factor1*(sin(factor) + sin(2.*factor));

  return partial_eff_lo;
}


/*--------------------------------------------------------------------------------------------------*/

/** IWASAKI **/

/* these are all implementations of analytical expressions */
/* note that we ONLY know analytically a solution for the EOM for the lattice background field
   of the IWASAKI acion in the case we choose Option_B of the paper: hep-lat/9808007, 
   which in this case, actually coincides with the expression obtained by the alpha_collab. for 
   the plaquette Wilson action with sf */
/* Therefore, whenever we are not considering this case B,
   the minimal lattice action (Iwasaki) configuration V must be found only numerically */


/* expression taken from the CP-PACS paper == K (eq. II.14) */
double partial_lattice_lo_effective_iwasaki_action_sf_k(int t, double beta, double c0, double c1, double eta) {

  double pi;
  double factor1, factor;
  double partial_eff_lo;

  pi = acos(-1.);

  factor1 = 12.*(double)LX*(double)LX;

  factor = (1./((double)LX*(double)t))*(eta + pi/3.0);


  partial_eff_lo = factor1*(c0*(sin(factor) + sin(2.*factor)) + 4.*c1*(sin(2.*factor) + sin(4.*factor)));

  return partial_eff_lo;
}

/*-------------------------------------------------------------------------------------------------*/

/*** DEFINITION OF THE RUNNING COUPLING ***/


/* it calculates an su3 matrix "u" which is gonna be the eighth-generator of SU(3) "lambda_8" */
#define _su3_lambda_8(u)					\
  (u).c00 = 1.0;						\
  (u).c01 = 0.0;						\
  (u).c02 = 0.0;						\
  (u).c10 = 0.0;						\
  (u).c11 =-0.5;						\
  (u).c12 = 0.0;						\
  (u).c20 = 0.0;						\
  (u).c21 = 0.0;						\
  (u).c22 =-0.5;						\


/* it calculates an su3 matrix "u" which is gonna be the eighth-generator of SU(3) "lambda_8" times "i" */
#define _su3_i_times_lambda_8(u)				\
  (u).c00 = 1.0 * I;						\
  (u).c01 = 0.0;						\
  (u).c02 = 0.0;						\
  (u).c10 = 0.0;						\
  (u).c11 =-0.5 * I;						\
  (u).c12 = 0.0;						\
  (u).c20 = 0.0;						\
  (u).c21 = 0.0;						\
  (u).c22 =-0.5 * I;

/*------------------------------------------------------------------------------------------------------------*/


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
	  _su3_times_su3(pr1,pr11,i_lambda8);
	  
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
  
  ga = ga_int;
  
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
	  _su3_times_su3(pr1,pr_r11,i_lambda8);	    
	  
	  
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
	  _su3_times_su3(pr1,pr_r11,i_lambda8);
	  
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

    ga = 2.*c1_tss*(ga_0_tss + ga_tminus1_tss) + c1_tts*(ga_0_tts + ga_tminus2_tts);


  return ga;  
  
  
}

/*------------------------------------------------------------------------------------------------------------*/

/* the next function gives us the "derivative" of the "WILSON action" with SF b.c.
   with respect to the background field parameter "eta" */
/* WARNING: this function is only valid if we are considering U!=V */
double partial_wilson_action_sf_respect_to_eta(int t, double beta, double cs, double ct) {

  double partial_plaquette;
  double partial_wilson;

  partial_plaquette = partial_plaquette_sf_respect_to_eta(t, ct);
  
  partial_wilson = - (beta/(3.*(double)LX)) * partial_plaquette;
  
  return partial_wilson;
}

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
  
  
  partial_iwasaki = - (beta/(3.*(double)LX)) * ( partial_plaquette + partial_rectangle );
  
  return partial_iwasaki;
}

