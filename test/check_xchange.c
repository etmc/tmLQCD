/*******************************************************************************
 * $Id$
 *
 * File check_xchange.c
 *
 * Check of the exchange routines
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "geometry_eo.h"
#ifdef MPI
#include "mpi_init.h"
#endif
#include "xchange.h"

void set_spinor_field(int k, const double c);
void set_gauge_field(const double c);
void set_spinor_point(spinor * s, const double c);
void set_deri_point();
su3 set_su3(const double c);
int check_geometry();

int main(int argc,char *argv[])
{
  double * x;
  int i,ix, mu, x0, x1, x2, x3;

  mpi_init(argc, argv);

  printf("\n");
  printf("Test of the mpi exchange routines \n");
  printf("----------------------------------\n");
  printf("\n");
  printf("The lattice size is %d x %d^3 \n\n",(int)(T*g_nproc_t),(int)(L));
   
  geometry();

  ix = check_geometry();

  /* Check the field exchange */
  /* Set the whole field to -1 */
  set_spinor_field(0, -1.);

  /* Set the internal boundary to g_cart_id */
  for(x1 = 0; x1 < LX; x1++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	set_spinor_point(&spinor_field[0][ g_ipt[0][x1][x2][x3]/2 ], g_cart_id);
	set_spinor_point(&spinor_field[0][ g_ipt[T-1][x1][x2][x3]/2 ], g_cart_id);
      }
    }
  }

#ifdef PARALLELXT
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	set_spinor_point(&spinor_field[0][ g_ipt[x0][0][x2][x3]/2 ], g_cart_id);
	set_spinor_point(&spinor_field[0][ g_ipt[x0][LX-1][x2][x3]/2 ], g_cart_id);
      }
    }
  }
#endif
  xchange_field(0);

#if (defined PARALLELT || defined PARALLELXT)  
  x = (double*) &spinor_field[0][T*LX*L*L/2];
  for(i = 0; i < LX*L*L/2*24; i++, x++) {
    if((int)(*x) != g_nb_t_up) {
      printf("The exchange up of fields in time direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
      printf("Program aborted\n");
#ifdef MPI
       MPI_Finalize(); 
#endif
       exit(0); 
    }
  }

  x = (double*) &spinor_field[0][(T+1)*LX*L*L/2];
  for(i = 0; i < LX*L*L/2*24; i++, x++) {
    if((int)(*x) != g_nb_t_dn) {
      printf("The exchange down of fields in time direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize(); 
#endif
      exit(0); 
    }
  }
#endif

#ifdef PARALLELXT
  x = (double*) &spinor_field[0][(T+2)*LX*L*L/2];
  for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
    if((int)(*x) != g_nb_x_up) {
      printf("The exchange up of fields in x direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize(); 
#endif
      exit(0); 
    }
  }

  x = (double*) &spinor_field[0][(T+2)*LX*L*L/2+T*LY*LZ/2];
  for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
    if((int)(*x) != g_nb_x_dn) {
      printf("The exchange down of fields in x direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize(); 
#endif
      exit(0); 
    }
  }
#endif

  /* Check the gauge exchange */

  set_gauge_field(-1.);

  for(x1 = 0; x1 < LX; x1++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu=0;mu<4;mu++){
	  g_gauge_field[ g_ipt[0][x1][x2][x3] ][mu]=set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][x2][x3] ][mu]=set_su3((double)g_cart_id);
	}
      }
    }
  }

#ifdef PARALLELXT
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu=0;mu<4;mu++){
	  g_gauge_field[ g_ipt[x0][0][x2][x3] ][mu]=set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][x2][x3] ][mu]=set_su3((double)g_cart_id);
	}
      }
    }
  }
#endif

  xchange_gauge();

#if (defined PARALLELT || defined PARALLELXT)  
  x = (double*) &g_gauge_field[T*LX*L*L][0];
  for(i = 0; i < LX*LY*LZ*72; i++, x++) {
    if((int)(*x) != g_nb_t_up) {
      printf("The exchange up of gaugefields in time direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

  x = (double*) &g_gauge_field[(T+1)*LX*L*L][0];
  for(i = 0; i < LX*LZ*LY*72; i++, x++) {
    if((int)(*x) != g_nb_t_dn) {
      printf("The exchange down of gaugefields in time direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }
#endif

#ifdef PARALLELXT
  x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ][0];
  for(i = 0; i < T*LY*LZ*72; i++, x++) {
    if((int)(*x) != g_nb_x_up) {
      printf("The exchange up of gaugefields in x direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
      printf("%d %d %d\n", g_cart_id, i, (int)(*x));
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

  x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ+T*LY*LZ][0];
  for(i = 0; i < T*LY*LZ*72; i++, x++) {
    if((int)(*x) != g_nb_x_dn) {
      printf("The exchange down of gaugefields in x direction\n");
      printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
      printf("%d %d %d\n", g_cart_id, i, (int)(*x));
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

#endif


  /* Check the deri exchange */

  for(ix = 0; ix < VOLUME+RAND; ix++){
    for(mu=0; mu<4; mu++){
      ddummy[ix][mu].d1=0.;
      ddummy[ix][mu].d2=0.;
      ddummy[ix][mu].d3=0.;
      ddummy[ix][mu].d4=0.;
      ddummy[ix][mu].d5=0.;
      ddummy[ix][mu].d6=0.;
      ddummy[ix][mu].d7=0.;
      ddummy[ix][mu].d8=0.;
      df0[ix][mu].d1=0.;
      df0[ix][mu].d2=0.;
      df0[ix][mu].d3=0.;
      df0[ix][mu].d4=0.;
      df0[ix][mu].d5=0.;
      df0[ix][mu].d6=0.;
      df0[ix][mu].d7=0.;
      df0[ix][mu].d8=0.;      
    }
  }

  for(x1 = 0; x1 < LX; x1++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T+1][x1][x2][x3];
	df0[ix][mu].d1=(double)g_cart_id;
	df0[ix][mu].d2=(double)g_cart_id;
	df0[ix][mu].d3=(double)g_cart_id;
	df0[ix][mu].d4=(double)g_cart_id;
	df0[ix][mu].d5=(double)g_cart_id;
	df0[ix][mu].d6=(double)g_cart_id;
	df0[ix][mu].d7=(double)g_cart_id;
	df0[ix][mu].d8=(double)g_cart_id;
      }
    }
  }
#ifdef PARALLELXT
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][LX+11][x2][x3];
	df0[ix][mu].d1=(double)g_cart_id;
	df0[ix][mu].d2=(double)g_cart_id;
	df0[ix][mu].d3=(double)g_cart_id;
	df0[ix][mu].d4=(double)g_cart_id;
	df0[ix][mu].d5=(double)g_cart_id;
	df0[ix][mu].d6=(double)g_cart_id;
	df0[ix][mu].d7=(double)g_cart_id;
	df0[ix][mu].d8=(double)g_cart_id;
      }
    }
  }
#endif

  xchange_deri();

  /* Check is missing! */

  printf("The exchange routines are working correctly\n");
  printf("\n");

#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}

void set_spinor_point(spinor * s, const double c){
  (*s).c1.c1.re=c;
  (*s).c1.c1.im=c;
  (*s).c1.c2.re=c;
  (*s).c1.c2.im=c;
  (*s).c1.c3.re=c;
  (*s).c1.c3.im=c;
  (*s).c2.c1.re=c;
  (*s).c2.c1.im=c;
  (*s).c2.c2.re=c;
  (*s).c2.c2.im=c;
  (*s).c2.c3.re=c;
  (*s).c2.c3.im=c;
  (*s).c3.c1.re=c;
  (*s).c3.c1.im=c;
  (*s).c3.c2.re=c;
  (*s).c3.c2.im=c;
  (*s).c3.c3.re=c;
  (*s).c3.c3.im=c;
  (*s).c4.c1.re=c;
  (*s).c4.c1.im=c;
  (*s).c4.c2.re=c;
  (*s).c4.c2.im=c;
  (*s).c4.c3.re=c;
  (*s).c4.c3.im=c;  
}

void set_spinor_field(int k, const double c) {

  int ix;
  spinor *s;
  for (ix=0;ix<VOLUME/2;ix++) {
    s=&spinor_field[k][ix];
    (*s).c1.c1.re=c;
    (*s).c1.c1.im=c;
    (*s).c1.c2.re=c;
    (*s).c1.c2.im=c;
    (*s).c1.c3.re=c;
    (*s).c1.c3.im=c;
    (*s).c2.c1.re=c;
    (*s).c2.c1.im=c;
    (*s).c2.c2.re=c;
    (*s).c2.c2.im=c;
    (*s).c2.c3.re=c;
    (*s).c2.c3.im=c;
    (*s).c3.c1.re=c;
    (*s).c3.c1.im=c;
    (*s).c3.c2.re=c;
    (*s).c3.c2.im=c;
    (*s).c3.c3.re=c;
    (*s).c3.c3.im=c;
    (*s).c4.c1.re=c;
    (*s).c4.c1.im=c;
    (*s).c4.c2.re=c;
    (*s).c4.c2.im=c;
    (*s).c4.c3.re=c;
    (*s).c4.c3.im=c;
  }
 for (ix=VOLUME/2;ix<VOLUMEPLUSRAND/2;ix++) {
    s=&spinor_field[k][ix];
    (*s).c1.c1.re=0;
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

su3 set_su3(const double c)
{
   su3 u;

   u.c11.re=c;
   u.c11.im=c;
   u.c12.re=c;
   u.c12.im=c;
   u.c13.re=c;
   u.c13.im=c;

   u.c21.re=c;
   u.c21.im=c;
   u.c22.re=c;
   u.c22.im=c;
   u.c23.re=c;
   u.c23.im=c;

   u.c31.re=c;
   u.c31.im=c;
   u.c32.re=c;
   u.c32.im=c;
   u.c33.re=c;
   u.c33.im=c;

   return(u);
}

void set_gauge_field(const double c) {
  int ix,mu;
  
  for (ix=0;ix<VOLUME;ix++) {
    for (mu=0;mu<4;mu++){
      g_gauge_field[ix][mu]=set_su3(c);
    }
  }
}

