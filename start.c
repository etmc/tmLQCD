/***********************************************************************
 *
 * Copyright (C) 2000 Martin Luescher
 *               2002 Martin Hasenbusch, Ines Wetzorke
 *               2003-2008 Carsten Urbach, Remi Baron
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
 * File start.c
 *
 * Collection of useful programs that mainly serve to initialize the fields
 *
 * The externally accessible functions are
 *
 *   su3_vector random_su3_vector(void)
 *     Returns a uniformly distributed random SU(3) vector with norm 1
 *
 *   spinor random_spinor(void)
 *     Returns a random spinor with norm 1
 *
 * M.Hasenbusch:
 *   void random_spinor_field(int k)
 *     Initializes the spinor field psi[k] to a Gaussian random field
 *
 * M.Hasenbusch:
 *   void zero_spinor_field(spinor * const k, const int V)
 *     Initializes the spinor field psi[k] to  zero
 *
 *   su3 random_su3(void)
 *     Returns a uniformly distributed random SU(3) matrix
 *
 *   void unit_g_gauge_field(void)
 *     Sets the gauge field variables to unity
 *
 *   void random_gauge_field(void)
 *     Initializes the gauge field to a random configuration
 *
 * Version: 1.0
 * Author: Martin Luescher <luscher@mail.desy.de>
 * Date: 24.10.2000
 *
 * Added the function
 *   void source_spinor_field_point_from_file(spinor * const P, spinor * const Q, int is, int ic, int source_indx)
 *   which uses the new input parameter SourceLocation in the input parameter files
 *   to place the source at the desired point
 *
 *   Author: Remi Baron <baron@th.u-psud.fr> April 2007
 *******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "read_input.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "ranlxs.h"
#include "start.h"

static void gauss_vector(double v[],int n)
{
   int k;
   double r[2];
   double x1,x2,rho,y1,y2;


   for (k=0;;k+=2)
   {
      ranlxd(r,2);
      x1=r[0];
      x2=r[1];

      rho = -log(1.0 - x1);
      rho = sqrt(rho);
      x2 *= 6.2831853071796;
      y1 = rho * sin(x2);
      y2 = rho * cos(x2);

      if (n > k)
         v[k] = y1;
      if (n > (k+1))
         v[k + 1] = y2;
      if (n <= (k + 2))
         return;
   }
}

/* produce a double array of z2 noise of length N */
static void z2_vector(double *v, const int N) {
  ranlxd(v,N);
  for (int i = 0; i < N; ++i) {
    if(v[i] < 0.5)
      v[i]=1/sqrt(2);
    else
      v[i]=-1/sqrt(2);
  }
  return;
}

static su3 unit_su3(void)
{
   su3 u = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
   return u;
}

su3_vector unit_su3_vector()
{
  su3_vector s = {1.0, 1.0, 1.0};
  return s;
}

/* produce a su3 vector with components distributed according to rn_type and unit norm */
static void random_su3_vector( su3_vector * const s, const enum RN_TYPE rn_type )
{
   int i;
   double v[6],norm,fact;

   void (*random_vector)(double*,int) = NULL;

   _rn_switch(rn_type,random_vector)

   while (1)
   {
      random_vector(v,6);
      norm=0.0;

      for (i = 0; i < 6; ++i)
         norm += v[i] * v[i];

      norm = sqrt(norm);

      if (1.0 != (1.0 + norm))
         break;
   }

   fact = 1.0 / norm;
   s->c0 = fact * (v[0] + I * v[1]);
   s->c1 = fact * (v[2] + I * v[3]);
   s->c2 = fact * (v[4] + I * v[5]);

   return;
}

static void random_spinor(spinor * const s, const enum RN_TYPE rn_type) {
   random_su3_vector(&s->s0, rn_type);
   random_su3_vector(&s->s1, rn_type);
   random_su3_vector(&s->s2, rn_type);
   random_su3_vector(&s->s3, rn_type);

   _vector_mul(s->s0, 0.5, s->s0);
   _vector_mul(s->s1, 0.5, s->s1);
   _vector_mul(s->s2, 0.5, s->s2);
   _vector_mul(s->s3, 0.5, s->s3);
   return;
}

spinor unit_spinor()
{
  spinor s;

  s.s0 = unit_su3_vector();
  s.s1 = unit_su3_vector();
  s.s2 = unit_su3_vector();
  s.s3 = unit_su3_vector();

  return(s);
}

void unit_spinor_field(const int k)
{
  int i=0;
  spinor *s;

  s = &g_spinor_field[k][0];
  for(i = 0; i < VOLUME/2; i++, s++) {
    *s = unit_spinor();
  }
}

/* Function provides a spinor field of length VOLUME with
   distributions given by rn_type as defined in start.h */
void random_spinor_field_lexic(spinor * const k, const int repro, const enum RN_TYPE rn_type) {
  int x, y, z, t, X, Y, Z, tt, id=0;

  void (*random_vector)(double*,int) = NULL;

  _rn_switch(rn_type,random_vector)

#ifdef MPI
  int rlxd_state[105];
  int rlxd_state_backup[105];
#endif
  int coords[4];
  spinor *s;
  double v[24];

  if(repro) {
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_get(rlxd_state_backup);
    } else if(g_proc_id == 0) {
      rlxd_get(rlxd_state);
    }
    MPI_Bcast(rlxd_state, 105, MPI_INT, 0, MPI_COMM_WORLD);
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state);
    }
#endif
    for(t = 0; t < g_nproc_t*T; t++) {
      tt = t - g_proc_coords[0]*T;
      coords[0] = t / T;
      for(x = 0; x < g_nproc_x*LX; x++) {
	X = x - g_proc_coords[1]*LX; 
	coords[1] = x / LX;
	for(y = 0; y < g_nproc_y*LY; y++) {
	  Y = y - g_proc_coords[2]*LY;
	  coords[2] = y / LY;
	  for(z = 0; z < g_nproc_z*LZ; z++) {
	    Z = z - g_proc_coords[3]*LZ;
	    coords[3] = z / LZ;
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    if(g_cart_id == id) {
	      random_vector(v, 24);
	      s = k + g_ipt[tt][X][Y][Z];
	      memcpy(s, v, 24*sizeof(double));
	    } else {
	      ranlxd(v,24);
	    }
	  }
	}
      }
    }
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state_backup);
    }
#endif
  }
  else {
    for(x = 0; x < VOLUME; x++) {
      random_vector(v, 24);
      s = k + x;
      memcpy(s, v, 24*sizeof(double));
    }
  }
  return;
}

/* Function provides a spinor field of length VOLUME/2 for even odd preconditioning 
   with distributions given by rn_type as defined in start.h */

void random_spinor_field_eo(spinor * const k, const int repro, const enum RN_TYPE rn_type ) {
  int x, X, y, Y, z, Z, t, t0, id = 0;

  void (*random_vector)(double*,int) = NULL;

  _rn_switch(rn_type,random_vector)

#ifdef MPI
  int rlxd_state[105];
  int rlxd_state_backup[105];
#endif
  int coords[4];
  spinor *s;
  double v[24];

  if(repro) {
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_get(rlxd_state_backup);
    } else if(g_proc_id == 0) {
      rlxd_get(rlxd_state);
    }
    MPI_Bcast(rlxd_state, 105, MPI_INT, 0, MPI_COMM_WORLD);
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state);
    }
#endif
    for(t0 = 0; t0 < g_nproc_t*T; t0++) {
      coords[0] = t0 / T;
      t = t0 - T*g_proc_coords[0];
      for(x = 0; x < g_nproc_x*LX; x++) {
	coords[1] = x / LX;
	X = x - g_proc_coords[1]*LX;
	for(y = 0; y < g_nproc_y*LY; y++) {
	  coords[2] = y / LY;
	  Y = y - g_proc_coords[2]*LY;
	  for(z = 0; z < g_nproc_z*LZ; z++) {
	    coords[3] = z / LZ;
	    Z = z - g_proc_coords[3]*LZ;
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    if((t0+x+y+z)%2 == 0) {
	      random_vector(v, 24);
	      if(g_cart_id == id) {
		s = k + g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
		memcpy(s, v, 24*sizeof(double));
	      }
	    }
	  }
	}
      }
    }
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state_backup);
    }
#endif
  }
  else {
    for (x = 0; x < VOLUME/2; x++) {
      s = k + x;
      random_vector(v, 24);
      memcpy(s, v, 24*sizeof(double));
    }
  }
  return;
}

/* Function provides a zero spinor field of length N */
void zero_spinor_field(spinor * const k, const int N)
{
  memset(k, 0, sizeof(spinor) * N);
}

/* Function provides a constant spinor field of length N */
void constant_spinor_field(spinor * const k, const int p, const int N)
{
  int ix;
  spinor *s;
  double * tmp;
  s = k;
  for (ix = 0; ix < N; ix++)
  {
    memset(s, 0, sizeof(spinor));
    tmp = (double*) s;
    tmp[2*p] = 1.;
    s++;
  }
  return;
}

/* a random su3 matrix. two unit norm su3 vectors with components drawn from a uniform ditribution
   are used to construct an orthogonal third vector. The three vectors then make up the rows
   of the matrix */
   
void random_su3(su3 * const u) {
   double norm,fact;
   _Complex double z;
   su3_vector z1,z2,z3;

   random_su3_vector(&z1,RN_UNIF);
   for (;;)
   {
      random_su3_vector(&z2,RN_UNIF);

      z = conj(z1.c0) * z2.c0 + conj(z1.c1) * z2.c1 + conj(z1.c2) * z2.c2;

      _vector_project(z2,z,z1);

      norm=sqrt(_vector_norm_square(z2));

      if (1.0 != (1.0 + norm))
         break;
   }

   fact = 1.0 / norm;
   _vector_mul(z2, fact, z2);

   z3.c0 = conj((z1.c1 * z2.c2) - (z1.c2 * z2.c1));
   z3.c1 = conj((z1.c2 * z2.c0) - (z1.c0 * z2.c2));
   z3.c2 = conj((z1.c0 * z2.c1) - (z1.c1 * z2.c0));

   u->c00 = z1.c0;
   u->c01 = z1.c1;
   u->c02 = z1.c2;

   u->c10 = z2.c0;
   u->c11 = z2.c1;
   u->c12 = z2.c2;

   u->c20 = z3.c0;
   u->c21 = z3.c1;
   u->c22 = z3.c2;
   return;
}


void unit_g_gauge_field(void)
{
  int ix,mu;

  for (ix=0;ix<VOLUME;ix++) {
    for (mu=0;mu<4;mu++) {
      g_gauge_field[ix][mu]=unit_su3();
    }
  }
  g_update_gauge_copy = 1;
  return;
}


void random_gauge_field(const int repro, su3 ** const gf) {

  int ix, mu, t0, t, x, X, y, Y, z, Z;
  int id = 0; /* May not be initialized for scalar builds! */
  int coords[4];
  su3 ALIGN tmp;
#ifdef MPI
  int rlxd_state[105];
  int rlxd_state_backup[105];
#endif

  if(repro) {
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_get(rlxd_state_backup);
    } else if(g_proc_id == 0) {
      rlxd_get(rlxd_state);
    }
    MPI_Bcast(rlxd_state, 105, MPI_INT, 0, MPI_COMM_WORLD);
    rlxd_reset(rlxd_state);
#endif
    for(t0 = 0; t0 < g_nproc_t*T; t0++) {
      t = t0 - T*g_proc_coords[0];
      coords[0] = t0 / T;
      for(x = 0; x < g_nproc_x*LX; x++) {
	X = x - g_proc_coords[1]*LX;
	coords[1] = x / LX;
	for(y = 0; y < g_nproc_y*LY; y++) {
	  Y = y - g_proc_coords[2]*LY;
	  coords[2] = y / LY;
	  for(z = 0; z < g_nproc_z*LZ; z++) {
	    Z = z - g_proc_coords[3]*LZ;
	    coords[3] = z / LZ;
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    for(mu = 0; mu < 4; mu++) {
	      if(g_cart_id == id) {
		ix = g_ipt[t][X][Y][Z];
		random_su3(&gf[ix][mu]);
	      }
	      else {
		random_su3(&tmp);
	      }
	    }
	  }
	}
      }
    }
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state_backup);
    }
#endif
  }
  else {
    for (ix = 0; ix < VOLUME; ix++) {
      for (mu = 0; mu < 4; mu++) {
	random_su3(&gf[ix][mu]);
      }
    }
  }

  g_update_gauge_copy = 1;
  return;
}

/* writes gaussian distributed random momenta of length VOLUME into momenta array
   and returns their energy contribution */
double random_su3adj_field(const int repro, su3adj ** const momenta) {
  su3adj *xm;
  int i, mu, t0, x, y, z, X, Y, Z, t, id = 0;
  int coords[4];
#ifdef MPI
  int k;
  int rlxd_state[105];
  int rlxd_state_backup[105];
#endif
  double ALIGN yy[8];
  double ALIGN tt, tr, ts, kc = 0., ks = 0., sum;
  
  if(repro) {
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_get(rlxd_state_backup);
    } else if(g_proc_id == 0) {
      rlxd_get(rlxd_state);
    }
    MPI_Bcast(rlxd_state, 105, MPI_INT, 0, MPI_COMM_WORLD);
    rlxd_reset(rlxd_state);
#endif
    for(t0 = 0; t0 < g_nproc_t*T; t0++) {
      t = t0 - T*g_proc_coords[0];
      coords[0] = t0 / T;
      for(x = 0; x < g_nproc_x*LX; x++) {
	X = x - g_proc_coords[1]*LX;
	coords[1] = x / LX;
	for(y = 0; y < g_nproc_y*LY; y++) {
	  Y = y - g_proc_coords[2]*LY;
	  coords[2] = y / LY;
	  for(z = 0; z < g_nproc_z*LZ; z++) {
	    Z = z - g_proc_coords[3]*LZ;
	    coords[3] = z / LZ;
#ifdef MPI
	    MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	    if(g_cart_id == id) i = g_ipt[t][X][Y][Z];
	    for(mu = 0; mu < 4; mu++) {
	      gauss_vector(yy,8);
	      if(g_cart_id == id) {
		sum = 0.;
		xm = &momenta[i][mu];
		(*xm).d1 = 1.4142135623731*yy[0];
		(*xm).d2 = 1.4142135623731*yy[1];
		sum += (*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
		(*xm).d3 = 1.4142135623731*yy[2];
		(*xm).d4 = 1.4142135623731*yy[3];
		sum += (*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
		(*xm).d5 = 1.4142135623731*yy[4];
		(*xm).d6 = 1.4142135623731*yy[5];
		sum += (*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	  (*xm).d7 = 1.4142135623731*yy[6];
	  (*xm).d8 = 1.4142135623731*yy[7];
	  sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	  tr=sum+kc;
		tr = sum+kc;
		ts = tr+ks;
		tt = ts-ks;
		ks = ts;
		kc = tr-tt;
	      }
	    }
	  }
	}
      }
    }
    kc=0.5*(ks+kc);
#ifdef MPI
    if(g_proc_id != 0) {
      rlxd_reset(rlxd_state_backup);
    }
#endif
  }
  else {
    for(i = 0; i < VOLUME; i++) { 
      for(mu = 0; mu < 4; mu++) {
	sum=0.;
	xm=&momenta[i][mu];
	gauss_vector(yy,8);
	(*xm).d1=1.4142135623731*yy[0];
	(*xm).d2=1.4142135623731*yy[1];
	sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	(*xm).d3=1.4142135623731*yy[2];
	(*xm).d4=1.4142135623731*yy[3];
	sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	(*xm).d5=1.4142135623731*yy[4];
	(*xm).d6=1.4142135623731*yy[5];
	sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	(*xm).d7=1.4142135623731*yy[6];
	(*xm).d8=1.4142135623731*yy[7];
	sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	tr=sum+kc;
	ts=tr+ks;
	tt=ts-ks;
	ks=ts;
	kc=tr-tt;
      }
    }
    kc=0.5*(ks+kc);
  }
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#endif
  return kc;
}

void set_spinor_point(spinor * s, const double c)
{
  s->s0.c0 = c * (1 + I);
  s->s0.c1 = c * (1 + I);
  s->s0.c2 = c * (1 + I);
  s->s1.c0 = c * (1 + I);
  s->s1.c1 = c * (1 + I);
  s->s1.c2 = c * (1 + I);
  s->s2.c0 = c * (1 + I);
  s->s2.c1 = c * (1 + I);
  s->s2.c2 = c * (1 + I);
  s->s3.c0 = c * (1 + I);
  s->s3.c1 = c * (1 + I);
  s->s3.c2 = c * (1 + I);
}

void set_spinor_field(int k, const double c)
{
  int ix;
  spinor *s;
  for (ix=0;ix<VOLUME/2;ix++)
  {
    s=&g_spinor_field[k][ix];
    s->s0.c0 = c * (1 + I);
    s->s0.c1 = c * (1 + I);
    s->s0.c2 = c * (1 + I);
    s->s1.c0 = c * (1 + I);
    s->s1.c1 = c * (1 + I);
    s->s1.c2 = c * (1 + I);
    s->s2.c0 = c * (1 + I);
    s->s2.c1 = c * (1 + I);
    s->s2.c2 = c * (1 + I);
    s->s3.c0 = c * (1 + I);
    s->s3.c1 = c * (1 + I);
    s->s3.c2 = c * (1 + I);
 }
 for (ix=VOLUME/2;ix<VOLUMEPLUSRAND/2;ix++)
 {
    s=&g_spinor_field[k][ix];
    s->s0.c0 = 0.;
    s->s0.c1 = 0.;
    s->s0.c2 = 0.;
    s->s1.c0 = 0.;
    s->s1.c1 = 0.;
    s->s1.c2 = 0.;
    s->s2.c0 = 0.;
    s->s2.c1 = 0.;
    s->s2.c2 = 0.;
    s->s3.c0 = 0.;
    s->s3.c1 = 0.;
    s->s3.c2 = 0.;
  }
}

su3 set_su3(const double c)
{
   su3 u;

   u.c00 = c * (1 + I);
   u.c01 = c * (1 + I);
   u.c02 = c * (1 + I);

   u.c10 = c * (1 + I);
   u.c11 = c * (1 + I);
   u.c12 = c * (1 + I);

   u.c20 = c * (1 + I);
   u.c21 = c * (1 + I);
   u.c22 = c * (1 + I);

   return(u);
}

void set_gauge_field(const double c)
{
  int ix,mu;

  for (ix=0;ix<VOLUMEPLUSRAND + g_dbw2rand;ix++) {
    for (mu=0;mu<4;mu++){
      g_gauge_field[ix][mu]=set_su3(c);
    }
  }
  g_update_gauge_copy = 1;
  return;
}


void source_spinor_field(spinor * const P, spinor * const Q, int is, int ic) {

  spinor * s;

  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  if (g_proc_coords[0] == 0 && g_proc_coords[1] == 0
      && g_proc_coords[2] == 0 && g_proc_coords[3] == 0) {

    s = P;

    /* put source to 1.0 */
    if (is==0){
      if      (ic==0) s->s0.c0 = 1.0;
      else if (ic==1) s->s0.c1 = 1.0;
      else if (ic==2) s->s0.c2 = 1.0;
    }
    else if (is==1){
      if      (ic==0) s->s1.c0 = 1.0;
      else if (ic==1) s->s1.c1 = 1.0;
      else if (ic==2) s->s1.c2 = 1.0;
    }
    else if (is==2){
      if      (ic==0) s->s2.c0 = 1.0;
      else if (ic==1) s->s2.c1 = 1.0;
      else if (ic==2) s->s2.c2 = 1.0;
    }
    else if (is==3){
      if      (ic==0) s->s3.c0 = 1.0;
      else if (ic==1) s->s3.c1 = 1.0;
      else if (ic==2) s->s3.c2 = 1.0;
    }
  }
}

void source_spinor_field_point_from_file(spinor * const P, spinor * const Q, int is, int ic, int source_indx)
{
  int tmp;
  int source_coord[4],source_pe_coord[4],source_loc_coord[4];
  int source_pe_indx,source_loc_indx;
  spinor * s;

  /* set fields to zero */
  zero_spinor_field(P,VOLUME/2);
  zero_spinor_field(Q,VOLUME/2);

  /* Check if source_indx is valid */
  if((source_indx < 0) || (source_indx >= (g_nproc_t*g_nproc_x*g_nproc_y*g_nproc_z*T*LX*LY*LZ)))
  {
    printf("Error in the input parameter file, SourceLocation must be in [0,VOLUME-1]! Exiting...!\n");
    exit(1);
  }

  /* translate it into global coordinate */
  /* For a T*L^3 lattice then  L = g_nproc_z * LZ = g_nproc_y * LY = g_nproc_x * LX    */
  source_coord[3]=source_indx % (g_nproc_z * LZ);
  tmp = source_indx / (g_nproc_z * LZ);
  source_coord[2]=tmp % (g_nproc_y * LY);
  tmp = tmp / (g_nproc_y * LY);
  source_coord[1]=tmp % (g_nproc_x * LX);
  tmp = tmp / (g_nproc_x * LX);
  source_coord[0]=tmp;

  if(3*is+ic == index_start && g_proc_id == g_stdio_proc)
    printf("# The source site number is %i which corresponds to (t,x,y,z) = (%i,%i,%i,%i)\n",source_indx,source_coord[0],source_coord[1],source_coord[2],source_coord[3]);

  /* compute the coordinates and the index of the node*/
  /* be careful!!! nodes indices have different convention (see io.c)*/
  source_pe_coord[0] = source_coord[0]/T;
  source_pe_coord[1] = source_coord[1]/LX;
  source_pe_coord[2] = source_coord[2]/LY;
  source_pe_coord[3] = source_coord[3]/LZ;

#ifdef MPI
  MPI_Cart_rank(g_cart_grid, source_pe_coord, &source_pe_indx);
#else
  source_pe_indx=0;
#endif

  /* compute the local (inside the node) coordinates and index*/
  source_loc_coord[0] = source_coord[0] - source_pe_coord[0] * T;
  source_loc_coord[1] = source_coord[1] - source_pe_coord[1] * LX;
  source_loc_coord[2] = source_coord[2] - source_pe_coord[2] * LY;
  source_loc_coord[3] = source_coord[3] - source_pe_coord[3] * LZ;

  source_loc_indx=g_ipt[source_loc_coord[0]][source_loc_coord[1]][source_loc_coord[2]][source_loc_coord[3]];

  /* Essayer g_proc_id au lieu de g_cart_id */
  if(source_pe_indx == g_cart_id)
  {
    if(3*is + ic == index_start && g_debug_level > 1)
    {
      printf("g_cart_id =%i\n",g_cart_id);
      printf("source_loc_coord[0] = %i\n",source_loc_coord[0]);
      printf("source_loc_coord[1] = %i\n",source_loc_coord[1]);
      printf("source_loc_coord[2] = %i\n",source_loc_coord[2]);
      printf("source_loc_coord[3] = %i\n",source_loc_coord[3]);
      printf("source_loc_indx = %i\n",source_loc_indx);
    }
    /* Check which spinor field (even or odd) needs to be initialized */
    if(g_lexic2eo[source_loc_indx] < VOLUME/2)
      s = P + g_lexic2eo[source_loc_indx];
    else
      s = Q + g_lexic2eosub[source_loc_indx];

    /* put source to 1.0 */
    if (is==0){
      if      (ic==0) s->s0.c0 = 1.0;
      else if (ic==1) s->s0.c1 = 1.0;
      else if (ic==2) s->s0.c2 = 1.0;
    }
    else if (is==1){
      if      (ic==0) s->s1.c0 = 1.0;
      else if (ic==1) s->s1.c1 = 1.0;
      else if (ic==2) s->s1.c2 = 1.0;
    }
    else if (is==2){
      if      (ic==0) s->s2.c0 = 1.0;
      else if (ic==1) s->s2.c1 = 1.0;
      else if (ic==2) s->s2.c2 = 1.0;
    }
    else if (is==3){
      if      (ic==0) s->s3.c0 = 1.0;
      else if (ic==1) s->s3.c1 = 1.0;
      else if (ic==2) s->s3.c2 = 1.0;
    }
  }
}

void start_ranlux(int level, int seed)
{
   unsigned int max_seed,loc_seed;
   unsigned int step = g_proc_coords[0]*g_nproc_x*g_nproc_y*g_nproc_z +
     g_proc_coords[1]*g_nproc_y*g_nproc_z +
     g_proc_coords[2]*g_nproc_z + g_proc_coords[3];

   max_seed = 2147483647 / g_nproc;
   loc_seed = (seed + step*max_seed) % 2147483647;

   if(loc_seed == 0) loc_seed++;

   #ifdef MPI
   unsigned int * seeds = calloc(g_nproc,sizeof(unsigned int));
   if(seeds == NULL) fatal_error("Memory allocation for seeds buffer failed!","start_ranlux");  
   MPI_Gather(&loc_seed,1,MPI_UNSIGNED,seeds,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
   if(g_proc_id == 0) {
     for(int i = 0; i < g_nproc; ++i) {
       for(int j = i+1; j < g_nproc; ++j) {
         if( seeds[i] == seeds[j] ) {
           char error_message[100];
           snprintf(error_message,100,"Process %d and %d have the same seed. Aborting!",i,j);
           fatal_error(error_message,"start_ranlux");
         }
       }
     }
   }
   free(seeds);
   #endif 
 
   if(g_debug_level > 3) {
     printf("Local seed is %d  proc_id = %d\n", loc_seed, g_proc_id);
   }

   rlxs_init(level-1, loc_seed);
   rlxd_init(level, loc_seed);
}

void gen_test_spinor_field(spinor * const k, const int eoflag) {

  int ix,iy,effvol;
  spinor *s;
  double invind,invvol;

  if (eoflag==1) {
    effvol=VOLUME/2;
  }else{
    effvol=VOLUME;
  }

  invvol=1/(VOLUME*100);
  s = k;

  for(ix = 0; ix < effvol; ix++){
    if (eoflag==1) {
      iy=g_eo2lexic[ix];
    }else{
      iy=ix;
    }
    
    invind=(double)(((g_coord[iy][0]*g_nproc_x*LX + g_coord[iy][1])*g_nproc_y*LY + g_coord[iy][2])*g_nproc_z*LZ + g_coord[iy][3] + 1.0);
    invind=1.0/invind;
    s->s0.c0 = invind;
    s->s0.c1 = invind+invvol;
    s->s0.c2 = invind+invvol/2.0;
    s->s1.c0 = invind+invvol/3.0;
    s->s1.c1 = invind+invvol/4.0;
    s->s1.c2 = invind+invvol/5.0;
    s->s2.c0 = invind+invvol/6.0;
    s->s2.c1 = invind+invvol/7.0;
    s->s2.c2 = invind+invvol/8.0;
    s->s3.c0 = invind+invvol/9.0;
    s->s3.c1 = invind+invvol/10.0;
    s->s3.c2 = invind+invvol/11.0;
    s++;
  }

}

void write_test_spinor_field(spinor * const k, const int eoflag, char * postfix) {
  FILE * testout;
  char  filenames[50];
  int ix,iy,effvol;

  sprintf(filenames,"test_out.%.4d.",g_proc_id);
  strcat(filenames,postfix);
  testout=fopen(filenames,"w");

  if (eoflag==1) {
    effvol=VOLUME/2;
  }else{
    effvol=VOLUME;
  }

  for(ix = 0; ix < effvol; ix++){
    if (eoflag==1) {
      iy=g_eo2lexic[ix];
    }else{
      iy=ix;
    }
    fprintf(testout,"[%d,%d,%d,%d;0,0]:%e\n",g_coord[iy][0],g_coord[iy][1],g_coord[iy][2],g_coord[iy][3],creal((k[ix]).s0.c0));
  }
  fclose(testout);
}
