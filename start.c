
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
*   void zero_spinor_field(int k)
*     Initializes the spinor field psi[k] to  zero
*
*   su3 random_su3(void)
*     Returns a uniformly distributed random SU(3) matrix
*
*   void unit_gauge_field(void)
*     Sets the gauge field variables to unity
*
*   void random_gauge_field(void)
*     Initializes the gauge field to a random configuration
*
* Version: 1.0
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 24.10.2000
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3adj.h"
#include "ranlxs.h"
#include "global.h"

/* murks to make it double precision by M.Hasenbusch */ 
void gauss_vector(double v[],int n)
{
   int k;
   float r[4];
/*   double pi; */
   double x1,x2,rho,y1,y2;
   
 /*  pi=4.0*atan(1.0); */

   for (k=0;;k+=2)
   {
      ranlxs(r,4);
      x1=((double)r[0])+0.596046448e-7*r[1];
      x2=((double)r[2])+0.596046448e-7*r[3];

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

   u.c11.re=1.0;
   u.c11.im=0.0;
   u.c12.re=0.0;
   u.c12.im=0.0;
   u.c13.re=0.0;
   u.c13.im=0.0;

   u.c21.re=0.0;
   u.c21.im=0.0;
   u.c22.re=1.0;
   u.c22.im=0.0;
   u.c23.re=0.0;
   u.c23.im=0.0;

   u.c31.re=0.0;
   u.c31.im=0.0;
   u.c32.re=0.0;
   u.c32.im=0.0;
   u.c33.re=1.0;
   u.c33.im=0.0;

   return(u);
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
   s.c1.re=v[0]*fact;
   s.c1.im=v[1]*fact;   
   s.c2.re=v[2]*fact;
   s.c2.im=v[3]*fact;
   s.c3.re=v[4]*fact;
   s.c3.im=v[5]*fact;     

   return(s);
}


spinor random_spinor(void)
{
   spinor s;

   s.c1=random_su3_vector();
   s.c2=random_su3_vector();
   s.c3=random_su3_vector();
   s.c4=random_su3_vector();   

   _vector_mul(s.c1,0.5,s.c1);
   _vector_mul(s.c2,0.5,s.c2);
   _vector_mul(s.c3,0.5,s.c3);
   _vector_mul(s.c4,0.5,s.c4);   
   
   return(s);
}

/* Function provides a spinor field of length VOLUME/2 with
   Gaussian distribution */
void random_spinor_field(int k)
{
int j,ix;
int rlxs_state[25];
spinor *s;
double v[6];

if(myid==0)
  {
for (ix=0;ix<VOLUME/2;ix++)
  {
  s=&spinor_field[k][ix];
  gauss_vector(v,6);
  (*s).c1.c1.re=v[0];
  (*s).c1.c1.im=v[1];
  (*s).c1.c2.re=v[2];
  (*s).c1.c2.im=v[3];
  (*s).c1.c3.re=v[4];
  (*s).c1.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c2.c1.re=v[0];
  (*s).c2.c1.im=v[1];
  (*s).c2.c2.re=v[2];
  (*s).c2.c2.im=v[3];
  (*s).c2.c3.re=v[4];
  (*s).c2.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c3.c1.re=v[0];
  (*s).c3.c1.im=v[1];
  (*s).c3.c2.re=v[2];
  (*s).c3.c2.im=v[3];
  (*s).c3.c3.re=v[4];
  (*s).c3.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c4.c1.re=v[0];
  (*s).c4.c1.im=v[1];
  (*s).c4.c2.re=v[2];
  (*s).c4.c2.im=v[3];
  (*s).c4.c3.re=v[4];
  (*s).c4.c3.im=v[5];
  }
/* send the state for the random-number generator to 1 */
  rlxs_get(rlxs_state);
  MPI_Send(&rlxs_state[0], 25, MPI_INT, 1, 102, MPI_COMM_WORLD);
  }

if(myid!=0)
  {
  MPI_Recv(&rlxs_state[0], 25, MPI_INT, myid-1, 102, MPI_COMM_WORLD, &status);
  rlxs_reset(rlxs_state);
for (ix=0;ix<VOLUME/2;ix++)
  {
  s=&spinor_field[k][ix];
  gauss_vector(v,6);
  (*s).c1.c1.re=v[0];
  (*s).c1.c1.im=v[1];
  (*s).c1.c2.re=v[2];
  (*s).c1.c2.im=v[3];
  (*s).c1.c3.re=v[4];
  (*s).c1.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c2.c1.re=v[0];
  (*s).c2.c1.im=v[1];
  (*s).c2.c2.re=v[2];
  (*s).c2.c2.im=v[3];
  (*s).c2.c3.re=v[4];
  (*s).c2.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c3.c1.re=v[0];
  (*s).c3.c1.im=v[1];
  (*s).c3.c2.re=v[2];
  (*s).c3.c2.im=v[3];
  (*s).c3.c3.re=v[4];
  (*s).c3.c3.im=v[5];
  gauss_vector(v,6);
  (*s).c4.c1.re=v[0];
  (*s).c4.c1.im=v[1];
  (*s).c4.c2.re=v[2];
  (*s).c4.c2.im=v[3];
  (*s).c4.c3.re=v[4];
  (*s).c4.c3.im=v[5];
  }
/* send the state fo the random-number generator to k+1 */
  j=myid+1; if(j==numprocs) j=0;
  rlxs_get(rlxs_state);
  MPI_Send(&rlxs_state[0], 25, MPI_INT, j, 102, MPI_COMM_WORLD);
  }
if(myid==0)
  {
  MPI_Recv(&rlxs_state[0], 25, MPI_INT, numprocs-1, 102, MPI_COMM_WORLD, &status);
  rlxs_reset(rlxs_state);
  }
}

/* Function provides a spinor field of length VOLUME/2 with
   Gaussian distribution */
void zero_spinor_field(int k)
{
int ix;
spinor *s;
for (ix=0;ix<VOLUME/2;ix++)
  {
  s=&spinor_field[k][ix];
  (*s).c1.c1.re=0.;
  (*s).c1.c1.im=0.;
  (*s).c1.c2.re=0.;
  (*s).c1.c2.im=0.;
  (*s).c1.c3.re=0.;
  (*s).c1.c3.im=0.;
  (*s).c2.c1.re=0.;
  (*s).c2.c1.im=0.;
  (*s).c2.c2.re=0.;
  (*s).c2.c2.im=0.;
  (*s).c2.c3.re=0.;
  (*s).c2.c3.im=0.;
  (*s).c3.c1.re=0.;
  (*s).c3.c1.im=0.;
  (*s).c3.c2.re=0.;
  (*s).c3.c2.im=0.;
  (*s).c3.c3.re=0.;
  (*s).c3.c3.im=0.;
  (*s).c4.c1.re=0.;
  (*s).c4.c1.im=0.;
  (*s).c4.c2.re=0.;
  (*s).c4.c2.im=0.;
  (*s).c4.c3.re=0.;
  (*s).c4.c3.im=0.;
  }
}

su3 random_su3(void)
{
   double norm,fact;
   complex z;
   su3_vector z1,z2,z3;
   su3 u;

   z1=random_su3_vector();

   for (;;)
   {
      z2=random_su3_vector();

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
   
   z3.c1.re= (z1.c2.re*z2.c3.re-z1.c2.im*z2.c3.im)
            -(z1.c3.re*z2.c2.re-z1.c3.im*z2.c2.im);
   z3.c1.im=-(z1.c2.re*z2.c3.im+z1.c2.im*z2.c3.re)
            +(z1.c3.re*z2.c2.im+z1.c3.im*z2.c2.re);

   z3.c2.re= (z1.c3.re*z2.c1.re-z1.c3.im*z2.c1.im)
            -(z1.c1.re*z2.c3.re-z1.c1.im*z2.c3.im);
   z3.c2.im=-(z1.c3.re*z2.c1.im+z1.c3.im*z2.c1.re)
            +(z1.c1.re*z2.c3.im+z1.c1.im*z2.c3.re);

   z3.c3.re= (z1.c1.re*z2.c2.re-z1.c1.im*z2.c2.im)
            -(z1.c2.re*z2.c1.re-z1.c2.im*z2.c1.im);
   z3.c3.im=-(z1.c1.re*z2.c2.im+z1.c1.im*z2.c2.re)
            +(z1.c2.re*z2.c1.im+z1.c2.im*z2.c1.re);    

   u.c11=z1.c1;
   u.c12=z1.c2;
   u.c13=z1.c3;

   u.c21=z2.c1;
   u.c22=z2.c2;
   u.c23=z2.c3;

   u.c31=z3.c1;
   u.c32=z3.c2;
   u.c33=z3.c3;   

   return(u);
}


void unit_gauge_field(void)
{
   int ix,mu;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         gauge_field[ix][mu]=unit_su3();
      }
   }
}


void random_gauge_field(void)
{
   int ix,mu;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         gauge_field[ix][mu]=random_su3();
      }
   }
}

