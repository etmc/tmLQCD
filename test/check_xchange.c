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
#include "start.h"
#include "xchange.h"

void set_deri_point();
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

  for(x1 = 0; x1 < 0; x1++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T+1][x1][x2][x3];
	for(mu=0;mu<4;mu++){
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
  }
#ifdef PARALLELXT
  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_idn[ g_ipt[x0][0][x2][x3] ][1];
	for(mu=0;mu<4;mu++){
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
  }
#endif

  xchange_deri();

  for(x0 = 0; x0 < T; x0++) {
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][LX-1][x2][x3];
	for(mu=0;mu<4;mu++){
	  if(
	  df0[ix][mu].d1 != g_nb_x_up ||
	  df0[ix][mu].d2 != g_nb_x_up ||
	  df0[ix][mu].d3 != g_nb_x_up ||
	  df0[ix][mu].d4 != g_nb_x_up ||
	  df0[ix][mu].d5 != g_nb_x_up ||
	  df0[ix][mu].d6 != g_nb_x_up ||
	  df0[ix][mu].d7 != g_nb_x_up ||
	  df0[ix][mu].d8 != g_nb_x_up){
	    printf("Exchange of derivatives is working not correctly!\n");
	    printf("Aborting program!");
#ifdef MPI
	    MPI_Finalize();
#endif
	    exit(0);
	  }
	}
      }
    }
  }

  /* Check is missing! */

  printf("The exchange routines are working correctly\n");
  printf("\n");

#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}


