/* $Id$ */

/*******************************************************************************
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
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "read_input.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "start.h"

void gauss_vector(double v[],int n)
{
   int k;
   double r[2];
/*    float r[4]; */
/*   double pi; */
   double x1,x2,rho,y1,y2;

 /*  pi=4.0*atan(1.0); */

   for (k=0;;k+=2)
   {
      ranlxd(r,2);
      x1=r[0];
      x2=r[1];

      rho=-log(1.0-x1);
      rho=sqrt(rho);
/*      x2*=2.0*pi; */
      x2*=6.2831853071796;
      y1=rho*sin(x2);
      y2=rho*cos(x2);

      if (n>k)
         v[k]=y1;
      if (n>(k+1))
         v[k+1]=y2;
      if (n<=(k+2))
         return;
   }
}


static su3 unit_su3(void)
{
   su3 u;

   u.c00.re=1.0;
   u.c00.im=0.0;
   u.c01.re=0.0;
   u.c01.im=0.0;
   u.c02.re=0.0;
   u.c02.im=0.0;

   u.c10.re=0.0;
   u.c10.im=0.0;
   u.c11.re=1.0;
   u.c11.im=0.0;
   u.c12.re=0.0;
   u.c12.im=0.0;

   u.c20.re=0.0;
   u.c20.im=0.0;
   u.c21.re=0.0;
   u.c21.im=0.0;
   u.c22.re=1.0;
   u.c22.im=0.0;

   return(u);
}

su3_vector unit_su3_vector() {
  su3_vector s;

  s.c0.re = 1.;
  s.c0.im = 0.;
  s.c1.re = 1.;
  s.c1.im = 0.;
  s.c2.re = 1.;
  s.c2.im = 0.;

  return(s);
}


su3_vector random_su3_vector(void)
{
   int i;
   double v[6],norm,fact;
   su3_vector s;

   for (;;)
   {
      gauss_vector(v,6);
      norm=0.0;

      for (i=0;i<6;i++)
         norm+=v[i]*v[i];

      norm=sqrt(norm);

      if (1.0!=(1.0+norm))
         break;
   }

   fact=1.0/norm;
   s.c0.re=v[0]*fact;
   s.c0.im=v[1]*fact;
   s.c1.re=v[2]*fact;
   s.c1.im=v[3]*fact;
   s.c2.re=v[4]*fact;
   s.c2.im=v[5]*fact;

   return(s);
}

su3_vector unif_su3_vector(void)
{
   int i;
   double v[6],norm,fact;
   su3_vector s;

   for (;;)
   {
      ranlxd(v,6);
      norm=0.0;

      for (i=0;i<6;i++){
	v[i] *= 6.2831853071796;
        norm+=v[i]*v[i];
     }

      norm=sqrt(norm);

      if (1.0!=(1.0+norm))
         break;
   }

   fact=1.0/norm;
   s.c0.re=v[0]*fact;
   s.c0.im=v[1]*fact;
   s.c1.re=v[2]*fact;
   s.c1.im=v[3]*fact;
   s.c2.re=v[4]*fact;
   s.c2.im=v[5]*fact;

   return(s);
}


spinor random_spinor(void)
{
   spinor s;

   s.s0=random_su3_vector();
   s.s1=random_su3_vector();
   s.s2=random_su3_vector();
   s.s3=random_su3_vector();

   _vector_mul(s.s0,0.5,s.s0);
   _vector_mul(s.s1,0.5,s.s1);
   _vector_mul(s.s2,0.5,s.s2);
   _vector_mul(s.s3,0.5,s.s3);

   return(s);
}

spinor unit_spinor() {
  spinor s;

  s.s0 = unit_su3_vector();
  s.s1 = unit_su3_vector();
  s.s2 = unit_su3_vector();
  s.s3 = unit_su3_vector();

  return(s);
}

void unit_spinor_field(const int k) {
  int i=0;
  spinor *s;

  s = &g_spinor_field[k][0];
  for(i = 0; i < VOLUME/2; i++, s++) {
    (*s) = unit_spinor();
  }
}

/* Function provides a spinor field of length V with
   Gaussian distribution */
void random_spinor_field(spinor * const k, const int V, const int repro) {

  int ix;
  int rlxd_state[105];
  spinor *s;
  double v[6];
#ifdef MPI
  int j=0;
#endif

  if(g_proc_id==0 && repro == 1) {
    for (ix = 0; ix < V; ix++) {
      s = k + ix;
      gauss_vector(v,6);
      (*s).s0.c0.re=v[0];
      (*s).s0.c0.im=v[1];
      (*s).s0.c1.re=v[2];
      (*s).s0.c1.im=v[3];
      (*s).s0.c2.re=v[4];
      (*s).s0.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s1.c0.re=v[0];
      (*s).s1.c0.im=v[1];
      (*s).s1.c1.re=v[2];
      (*s).s1.c1.im=v[3];
      (*s).s1.c2.re=v[4];
      (*s).s1.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s2.c0.re=v[0];
      (*s).s2.c0.im=v[1];
      (*s).s2.c1.re=v[2];
      (*s).s2.c1.im=v[3];
      (*s).s2.c2.re=v[4];
      (*s).s2.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s3.c0.re=v[0];
      (*s).s3.c0.im=v[1];
      (*s).s3.c1.re=v[2];
      (*s).s3.c1.im=v[3];
      (*s).s3.c2.re=v[4];
      (*s).s3.c2.im=v[5];
    }
    /* send the state for the random-number generator to 1 */
    rlxd_get(rlxd_state);
#ifdef MPI
    MPI_Send(&rlxd_state[0], 105, MPI_INT, 1, 102, MPI_COMM_WORLD);
#endif
  }
#ifdef MPI
  if(g_proc_id != 0 && repro == 1) {
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_proc_id-1, 102, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
    for (ix=0;ix<V;ix++) {
      s = k + ix;
      gauss_vector(v,6);
      (*s).s0.c0.re=v[0];
      (*s).s0.c0.im=v[1];
      (*s).s0.c1.re=v[2];
      (*s).s0.c1.im=v[3];
      (*s).s0.c2.re=v[4];
      (*s).s0.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s1.c0.re=v[0];
      (*s).s1.c0.im=v[1];
      (*s).s1.c1.re=v[2];
      (*s).s1.c1.im=v[3];
      (*s).s1.c2.re=v[4];
      (*s).s1.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s2.c0.re=v[0];
      (*s).s2.c0.im=v[1];
      (*s).s2.c1.re=v[2];
      (*s).s2.c1.im=v[3];
      (*s).s2.c2.re=v[4];
      (*s).s2.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s3.c0.re=v[0];
      (*s).s3.c0.im=v[1];
      (*s).s3.c1.re=v[2];
      (*s).s3.c1.im=v[3];
      (*s).s3.c2.re=v[4];
      (*s).s3.c2.im=v[5];
    }
    /* send the state fo the random-number generator to k+1 */
    j=g_proc_id+1;
    if(j==g_nproc){
      j=0;
    }
    rlxd_get(rlxd_state);
    MPI_Send(&rlxd_state[0], 105, MPI_INT, j, 102, MPI_COMM_WORLD);
  }
  if(g_proc_id==0 && repro == 1) {
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 102, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
  }
#endif
  if(repro != 1) {
    for (ix = 0; ix < V; ix++) {
      s = k + ix;
      gauss_vector(v,6);
      (*s).s0.c0.re=v[0];
      (*s).s0.c0.im=v[1];
      (*s).s0.c1.re=v[2];
      (*s).s0.c1.im=v[3];
      (*s).s0.c2.re=v[4];
      (*s).s0.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s1.c0.re=v[0];
      (*s).s1.c0.im=v[1];
      (*s).s1.c1.re=v[2];
      (*s).s1.c1.im=v[3];
      (*s).s1.c2.re=v[4];
      (*s).s1.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s2.c0.re=v[0];
      (*s).s2.c0.im=v[1];
      (*s).s2.c1.re=v[2];
      (*s).s2.c1.im=v[3];
      (*s).s2.c2.re=v[4];
      (*s).s2.c2.im=v[5];
      gauss_vector(v,6);
      (*s).s3.c0.re=v[0];
      (*s).s3.c0.im=v[1];
      (*s).s3.c1.re=v[2];
      (*s).s3.c1.im=v[3];
      (*s).s3.c2.re=v[4];
      (*s).s3.c2.im=v[5];
    }
  }
}

/* Function provides a zero spinor field of length N with */
void zero_spinor_field(spinor * const k, const int N) {

  int ix;
  spinor *s;
  s = k;
  for (ix=0;ix<N;ix++) {
    (*s).s0.c0.re=0.;
    (*s).s0.c0.im=0.;
    (*s).s0.c1.re=0.;
    (*s).s0.c1.im=0.;
    (*s).s0.c2.re=0.;
    (*s).s0.c2.im=0.;
    (*s).s1.c0.re=0.;
    (*s).s1.c0.im=0.;
    (*s).s1.c1.re=0.;
    (*s).s1.c1.im=0.;
    (*s).s1.c2.re=0.;
    (*s).s1.c2.im=0.;
    (*s).s2.c0.re=0.;
    (*s).s2.c0.im=0.;
    (*s).s2.c1.re=0.;
    (*s).s2.c1.im=0.;
    (*s).s2.c2.re=0.;
    (*s).s2.c2.im=0.;
    (*s).s3.c0.re=0.;
    (*s).s3.c0.im=0.;
    (*s).s3.c1.re=0.;
    (*s).s3.c1.im=0.;
    (*s).s3.c2.re=0.;
    (*s).s3.c2.im=0.;
    s++;
  }
  return;
}

/* Function provides a constant spinor field of length N with */
void constant_spinor_field(spinor * const k, const int p, const int N) {

  int ix;
  spinor *s;
  double * tmp;
  s = k;
  for (ix = 0; ix < N; ix++) {
    (*s).s0.c0.re=0.;
    (*s).s0.c0.im=0.;
    (*s).s0.c1.re=0.;
    (*s).s0.c1.im=0.;
    (*s).s0.c2.re=0.;
    (*s).s0.c2.im=0.;
    (*s).s1.c0.re=0.;
    (*s).s1.c0.im=0.;
    (*s).s1.c1.re=0.;
    (*s).s1.c1.im=0.;
    (*s).s1.c2.re=0.;
    (*s).s1.c2.im=0.;
    (*s).s2.c0.re=0.;
    (*s).s2.c0.im=0.;
    (*s).s2.c1.re=0.;
    (*s).s2.c1.im=0.;
    (*s).s2.c2.re=0.;
    (*s).s2.c2.im=0.;
    (*s).s3.c0.re=0.;
    (*s).s3.c0.im=0.;
    (*s).s3.c1.re=0.;
    (*s).s3.c1.im=0.;
    (*s).s3.c2.re=0.;
    (*s).s3.c2.im=0.;
    tmp = (double*) s;
    tmp[2*p] = 1.;
    s++;
  }
  return;
}


su3 random_su3(void)
{
   double norm,fact;
   complex z;
   su3_vector z1,z2,z3;
   su3 u;

   /*
   z1=random_su3_vector();
   */
   z1=unif_su3_vector();


   for (;;)
   {
     /*
      z2=random_su3_vector();
     */
      z2=unif_su3_vector();

      z.re=_vector_prod_re(z1,z2);
      z.im=_vector_prod_im(z1,z2);

      _vector_project(z2,z,z1);

      norm=_vector_prod_re(z2,z2);
      norm=sqrt(norm);

      if (1.0!=(1.0+norm))
         break;
   }

   fact=1.0/norm;
   _vector_mul(z2,fact,z2);

   z3.c0.re= (z1.c1.re*z2.c2.re-z1.c1.im*z2.c2.im)
            -(z1.c2.re*z2.c1.re-z1.c2.im*z2.c1.im);
   z3.c0.im=-(z1.c1.re*z2.c2.im+z1.c1.im*z2.c2.re)
            +(z1.c2.re*z2.c1.im+z1.c2.im*z2.c1.re);

   z3.c1.re= (z1.c2.re*z2.c0.re-z1.c2.im*z2.c0.im)
            -(z1.c0.re*z2.c2.re-z1.c0.im*z2.c2.im);
   z3.c1.im=-(z1.c2.re*z2.c0.im+z1.c2.im*z2.c0.re)
            +(z1.c0.re*z2.c2.im+z1.c0.im*z2.c2.re);

   z3.c2.re= (z1.c0.re*z2.c1.re-z1.c0.im*z2.c1.im)
            -(z1.c1.re*z2.c0.re-z1.c1.im*z2.c0.im);
   z3.c2.im=-(z1.c0.re*z2.c1.im+z1.c0.im*z2.c1.re)
            +(z1.c1.re*z2.c0.im+z1.c1.im*z2.c0.re);

   u.c00=z1.c0;
   u.c01=z1.c1;
   u.c02=z1.c2;

   u.c10=z2.c0;
   u.c11=z2.c1;
   u.c12=z2.c2;

   u.c20=z3.c0;
   u.c21=z3.c1;
   u.c22=z3.c2;

   return(u);
}


void unit_g_gauge_field(void) {
  int ix,mu;

  for (ix=0;ix<VOLUME;ix++) {
    for (mu=0;mu<4;mu++) {
      g_gauge_field[ix][mu]=unit_su3();
    }
  }
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
  return;
}


void random_gauge_field() {

  int ix,mu;
#ifdef MPI
  int rlxd_state[105];
  int j=0;

  if(g_proc_id !=0) {
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_proc_id-1, 102, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
  }
#endif

  for (ix = 0; ix < VOLUME; ix++) {
    for (mu = 0; mu < 4; mu++) {
      g_gauge_field[ix][mu] = random_su3();
    }
  }

#ifdef MPI
  j = (g_proc_id + 1) % g_nproc;
  rlxd_get(rlxd_state);
  MPI_Send(&rlxd_state[0], 105, MPI_INT, j, 102, MPI_COMM_WORLD);

  if(g_proc_id == 0) {
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 102, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
  }
#endif
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
  return;
}

void set_spinor_point(spinor * s, const double c){
  (*s).s0.c0.re=c;
  (*s).s0.c0.im=c;
  (*s).s0.c1.re=c;
  (*s).s0.c1.im=c;
  (*s).s0.c2.re=c;
  (*s).s0.c2.im=c;
  (*s).s1.c0.re=c;
  (*s).s1.c0.im=c;
  (*s).s1.c1.re=c;
  (*s).s1.c1.im=c;
  (*s).s1.c2.re=c;
  (*s).s1.c2.im=c;
  (*s).s2.c0.re=c;
  (*s).s2.c0.im=c;
  (*s).s2.c1.re=c;
  (*s).s2.c1.im=c;
  (*s).s2.c2.re=c;
  (*s).s2.c2.im=c;
  (*s).s3.c0.re=c;
  (*s).s3.c0.im=c;
  (*s).s3.c1.re=c;
  (*s).s3.c1.im=c;
  (*s).s3.c2.re=c;
  (*s).s3.c2.im=c;
}

void set_spinor_field(int k, const double c) {

  int ix;
  spinor *s;
  for (ix=0;ix<VOLUME/2;ix++) {
    s=&g_spinor_field[k][ix];
    (*s).s0.c0.re=c;
    (*s).s0.c0.im=c;
    (*s).s0.c1.re=c;
    (*s).s0.c1.im=c;
    (*s).s0.c2.re=c;
    (*s).s0.c2.im=c;
    (*s).s1.c0.re=c;
    (*s).s1.c0.im=c;
    (*s).s1.c1.re=c;
    (*s).s1.c1.im=c;
    (*s).s1.c2.re=c;
    (*s).s1.c2.im=c;
    (*s).s2.c0.re=c;
    (*s).s2.c0.im=c;
    (*s).s2.c1.re=c;
    (*s).s2.c1.im=c;
    (*s).s2.c2.re=c;
    (*s).s2.c2.im=c;
    (*s).s3.c0.re=c;
    (*s).s3.c0.im=c;
    (*s).s3.c1.re=c;
    (*s).s3.c1.im=c;
    (*s).s3.c2.re=c;
    (*s).s3.c2.im=c;
  }
 for (ix=VOLUME/2;ix<VOLUMEPLUSRAND/2;ix++) {
    s=&g_spinor_field[k][ix];
    (*s).s0.c0.re=0;
    (*s).s0.c0.im=0.;
    (*s).s0.c1.re=0.;
    (*s).s0.c1.im=0.;
    (*s).s0.c2.re=0.;
    (*s).s0.c2.im=0.;
    (*s).s1.c0.re=0.;
    (*s).s1.c0.im=0.;
    (*s).s1.c1.re=0.;
    (*s).s1.c1.im=0.;
    (*s).s1.c2.re=0.;
    (*s).s1.c2.im=0.;
    (*s).s2.c0.re=0.;
    (*s).s2.c0.im=0.;
    (*s).s2.c1.re=0.;
    (*s).s2.c1.im=0.;
    (*s).s2.c2.re=0.;
    (*s).s2.c2.im=0.;
    (*s).s3.c0.re=0.;
    (*s).s3.c0.im=0.;
    (*s).s3.c1.re=0.;
    (*s).s3.c1.im=0.;
    (*s).s3.c2.re=0.;
    (*s).s3.c2.im=0.;
  }
}

su3 set_su3(const double c)
{
   su3 u;

   u.c00.re=c;
   u.c00.im=c;
   u.c01.re=c;
   u.c01.im=c;
   u.c02.re=c;
   u.c02.im=c;

   u.c10.re=c;
   u.c10.im=c;
   u.c11.re=c;
   u.c11.im=c;
   u.c12.re=c;
   u.c12.im=c;

   u.c20.re=c;
   u.c20.im=c;
   u.c21.re=c;
   u.c21.im=c;
   u.c22.re=c;
   u.c22.im=c;

   return(u);
}

void set_gauge_field(const double c) {
  int ix,mu;

  for (ix=0;ix<VOLUMEPLUSRAND + g_dbw2rand;ix++) {
    for (mu=0;mu<4;mu++){
      g_gauge_field[ix][mu]=set_su3(c);
    }
  }
  g_update_gauge_copy = 1;
  g_update_gauge_energy = 1;
  g_update_rectangle_energy = 1;
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
      if      (ic==0) (*s).s0.c0.re=1.0;
      else if (ic==1) (*s).s0.c1.re=1.0;
      else if (ic==2) (*s).s0.c2.re=1.0;
    }
    else if (is==1){
      if      (ic==0) (*s).s1.c0.re=1.0;
      else if (ic==1) (*s).s1.c1.re=1.0;
      else if (ic==2) (*s).s1.c2.re=1.0;
    }
    else if (is==2){
      if      (ic==0) (*s).s2.c0.re=1.0;
      else if (ic==1) (*s).s2.c1.re=1.0;
      else if (ic==2) (*s).s2.c2.re=1.0;
    }
    else if (is==3){
      if      (ic==0) (*s).s3.c0.re=1.0;
      else if (ic==1) (*s).s3.c1.re=1.0;
      else if (ic==2) (*s).s3.c2.re=1.0;
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
    printf("The source site number is %i which corresponds to (t,x,y,z)=(%i,%i,%i,%i)\n",source_indx,source_coord[0],source_coord[1],source_coord[2],source_coord[3]);

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
      if      (ic==0) (*s).s0.c0.re=1.0;
      else if (ic==1) (*s).s0.c1.re=1.0;
      else if (ic==2) (*s).s0.c2.re=1.0;
    }
    else if (is==1){
      if      (ic==0) (*s).s1.c0.re=1.0;
      else if (ic==1) (*s).s1.c1.re=1.0;
      else if (ic==2) (*s).s1.c2.re=1.0;
    }
    else if (is==2){
      if      (ic==0) (*s).s2.c0.re=1.0;
      else if (ic==1) (*s).s2.c1.re=1.0;
      else if (ic==2) (*s).s2.c2.re=1.0;
    }
    else if (is==3){
      if      (ic==0) (*s).s3.c0.re=1.0;
      else if (ic==1) (*s).s3.c1.re=1.0;
      else if (ic==2) (*s).s3.c2.re=1.0;
    }
  }
}
