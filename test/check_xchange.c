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

int main(int argc,char *argv[])
{
  double * x;
  int i,ix;

  mpi_init(argc, argv);

  printf("\n");
  printf("Test of the mpi exchange routines \n");
  printf("----------------------------------\n");
  printf("\n");
  printf("The lattice size is %d x %d^3 \n\n",(int)(T*g_nproc_t),(int)(L));
   
  geometry();

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
    }
  }

  /* Check the field exchange */
  set_spinor_field(0, (double)g_proc_id);

  xchange_field(0);

#if (defined PARALLELT || defined PARALLELXT)  
  x = (double*) &spinor_field[0][T*LX*L*L/2];
  for(i = 0; i < LX*L*L/2*12; i++, x++) {
    if((int)(*x) != g_nb_t_up) {
      printf("The exchange up of fields in time direction\n");
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_t_up);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

  x = (double*) &spinor_field[0][(T+1)*LX*L*L/2];
  for(i = 0; i < LX*L*L/2*12; i++, x++) {
    if((int)(*x) != g_nb_t_dn) {
      printf("The exchange down of fields in time direction\n");
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_t_dn);
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
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_x_up);
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
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_x_dn);
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }
#endif

  /* Check the gauge exchange */

  set_gauge_field((double)g_proc_id);

  xchange_gauge();

#if (defined PARALLELT || defined PARALLELXT)  
  x = (double*) &g_gauge_field[T*LX*L*L][0];
  for(i = 0; i < LX*LY*LZ*72; i++, x++) {
    if((int)(*x) != g_nb_t_up) {
      printf("The exchange up of gaugefields in time direction\n");
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_t_up);
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
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_t_dn);
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
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_x_up);
      printf("%d %d %d\n", g_proc_id, i, (int)(*x));
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
      printf("between %d and %d is not correct\n", g_proc_id, g_nb_x_dn);
      printf("%d %d %d\n", g_proc_id, i, (int)(*x));
      printf("Program aborted\n");
#ifdef MPI
      MPI_Finalize();
#endif
      exit(0);
    }
  }

#endif


  /* Check the deri exchange */

  xchange_deri();

  printf("The exchange routines are working correctly\n");
  printf("\n");

#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
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

